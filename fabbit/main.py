import argparse
from pathlib import Path
import os
from io import StringIO, BytesIO
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import pyrodigal
from Bio import SeqIO
import pyfastx
import logging
import traceback
import gzip
from fabbit.diamond_wrapper import DIAMOND
import pickle
import pandas as pd
from fabbit.mafft_wrapper import MAFFT
import numpy as np
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
import tempfile
import logging
import matplotlib.pyplot as plt
import time
import asyncio
from concurrent.futures import ThreadPoolExecutor
import gc
from tqdm.asyncio import tqdm as atqdm

__version__ = "1.0.0"  # Single source of truth for version

def setup_logging(verbose_level, output_dir):
    log_levels = {
        1: logging.ERROR,
        2: logging.WARNING,
        3: logging.INFO,
        4: logging.DEBUG
    }
    log_level = log_levels.get(verbose_level, logging.INFO)
    
    # Create logs directory if it doesn't exist
    log_dir = Path(output_dir) / 'logs'
    log_dir.mkdir(parents=True, exist_ok=True)
    
    # Set up logging to file
    log_file = log_dir / 'fabbit.log.txt'
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        filename=str(log_file),
        filemode='w'
    )
    
    # Optionally, add a console handler if you still want some minimal output to console
    console = logging.StreamHandler()
    console.setLevel(logging.ERROR)  # Only show errors in console
    formatter = logging.Formatter('%(levelname)s - %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

def parse_arguments():
    parser = argparse.ArgumentParser(description="FABBIT: FAst coregenome alignment Based on Bidirectional best hIT")
    parser.add_argument("-i", "--input", nargs='+', required=True, help="Path to the input FASTA file(s)")
    parser.add_argument("-o", "--output", required=True, help="Path to the output directory")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use")
    parser.add_argument("-r", "--reference", help="Path to the reference genome file (optional)")
    parser.add_argument("--sensitivity", default="fast", choices=["fast", "mid-sensitive", "sensitive", "more-sensitive", "very-sensitive", "ultra-sensitive"], help="DIAMOND sensitivity mode")
    parser.add_argument("--evalue", type=float, default=1e-6, help="Maximum e-value to report alignments")
    parser.add_argument("--query-cover", type=float, default=95, help="Minimum query cover percentage")
    parser.add_argument("--max-target-seqs", type=int, default=25, help="Maximum number of target sequences per query")
    parser.add_argument("--id", type=float, default=30, help="Minimum identity percentage")
    parser.add_argument("--core-threshold", type=float, default=95, help="Threshold for defining core genes")
    parser.add_argument("-v", "--verbose", type=int, choices=[1, 2, 3, 4], default=2,
                        help="Verbose level: 1=ERROR, 2=WARNING, 3=INFO, 4=DEBUG")
    # Add the version argument
    parser.add_argument("-V", "--version", action="version", version=f"%(prog)s {__version__}")
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    return args, output_dir

def calculate_shannon_entropy(align: MultipleSeqAlignment):
    # Transform the alignment into a compact numerical representation
    alignment_array = np.array([list(str(record.seq).lower()) for record in align], dtype='|S1').view(np.int8)
    
    # Define ASCII codes for nucleotides and the code for gap/other
    nucleotides = np.array([97, 99, 103, 116, 110], dtype=np.int8)  # a, c, g, t, n(others)
    
    # Mask for values different from a, c, g, t
    mask = ~np.isin(alignment_array, nucleotides[:4])
    alignment_array[mask] = 110  # assign 'n' to non-a,c,g,t
    
    # Count occurrences
    col_counts = np.array([np.sum(alignment_array == nt_code, axis=0) for nt_code in nucleotides], dtype=np.float64)
    
    # Normalize counts to get frequencies
    freqs = col_counts / np.sum(col_counts, axis=0)
    
    # Calculate Shannon entropy
    with np.errstate(divide='ignore', invalid='ignore'):  # Ignore log(0) errors
        entropies = -np.nansum(freqs[:4,:] * np.log2(freqs[:4,:]), axis=0)  # ignore 'n' in calculation
    return np.nanmean(entropies)  # Return the mean entropy of the columns

def create_output_directories(output_dir):
    subdirs = ['pyrodigal_orfs', 'core_genome_genes']
    for subdir in subdirs:
        path = os.path.join(output_dir, subdir)
        os.makedirs(path, exist_ok=True)
    logging.info(f"Created output directories: {', '.join(subdirs)}")
    return subdirs

def predict_orfs(filename, output_dir, meta_mode=False):
    genome_path = filename
    filebase = filename.stem
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    faa_path = Path(output_dir) / f"{filebase}.faa"
    fna_path = Path(output_dir) / f"{filebase}.fna"

    # If both files exist, return without doing anything
    if faa_path.exists() and fna_path.exists():
        logging.info(f"Skipping {filename}: Output files already exist")
        return

    logging.info(f"Processing {filename}")

    # Use StringIO objects to accumulate the outputs
    faa_output = StringIO()
    fna_output = StringIO()

    orf_finder = pyrodigal.GeneFinder(meta=meta_mode)

    # accumulated_seq to train prodigal
    accumulated_seq = BytesIO()

    # Determine if the file is gzipped
    is_gzipped = filename.suffix.lower() == '.gz'

    # Open the file with the appropriate method
    open_func = gzip.open if is_gzipped else open
    mode = 'rt' if is_gzipped else 'r'

    with open_func(genome_path, mode) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            accumulated_seq.write(bytes(record.seq))
            if len(accumulated_seq.getvalue()) >= 100000:
                break

    if not meta_mode:
        logging.info(f"Training ORF finder for {filename}")
        orf_finder.train(accumulated_seq.getvalue())
        meta_mode = True  # Prevent re-training if training once is enough

    # Process the records in the input file
    with open_func(genome_path, mode) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            genes = orf_finder.find_genes(bytes(record.seq))
            
            # Write translations (proteins) to StringIO object
            genes.write_translations(faa_output, sequence_id=record.id, width=80)
            
            # Write genes (nucleotide sequences) to StringIO object
            genes.write_genes(fna_output, sequence_id=record.id, width=80)
    
    # Write the accumulated outputs to disk from the StringIO objects
    with open(faa_path, 'w') as faa_file:
        faa_output.seek(0)
        faa_file.write(faa_output.read())
    
    with open(fna_path, 'w') as fna_file:
        fna_output.seek(0)
        fna_file.write(fna_output.read())

    # Close the StringIO objects
    faa_output.close()
    fna_output.close()
    pyfastx.Fasta(str(fna_path))
    logging.info(f"Finished processing {filename}")

def process_genomes(input_files, output_dir, max_workers=None):
    # Create the output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    logging.info(f"Processing {len(input_files)} genomes with {max_workers} workers")

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Create a tqdm iterable for a progress bar
        futures = {executor.submit(predict_orfs, filename, output_dir): filename for filename in input_files}
        
        for future in tqdm(as_completed(futures), total=len(input_files)):
            filename = futures[future]
            try:
                future.result()
            except Exception as e:
                logging.error(f"Exception occurred while processing {filename}: {traceback.format_exc()}")
                raise e

def count_fasta_headers(reference_file):
    """
    Counts the number of FASTA headers (lines starting with '>') in the reference file.
    """
    with open(reference_file, 'r') as file:
        content = file.read()
        header_count = content.count('>')
    return header_count

def parallel_bidirectional_best_hit(reference, queries, diamond_instance, output, core_threshold, **kwargs):
    """
    Perform Bidirectional Best Hit (BBH) analysis between a reference genome and multiple query genomes using parallel processing.

    Parameters:
    - reference (str or Path): Path to the reference genome FASTA file.
    - queries (list): List of paths to query genome FASTA files.
    - diamond_instance: Instance with a method 'bidirectional_best_hit' for performing BBH.
    - output (Path): Directory where output files will be saved.
    - core_threshold (float): Threshold (0-1) for defining core genes based on their presence across genomes.
    - **kwargs: Additional keyword arguments to pass to 'diamond_instance.bidirectional_best_hit'.

    Returns:
    - result_df (pd.DataFrame): Summary results for each genome.
    - coregenome_matrix (pd.DataFrame): Matrix representing the core genome across genomes.
    """
    # Ensure output directory exists
    output.mkdir(parents=True, exist_ok=True)

    # Extract 'threads' from kwargs or set default
    threads = kwargs.get('threads', os.cpu_count())
    # Remove 'threads' from kwargs to avoid passing it to 'bidirectional_best_hit'
    kwargs.pop('threads', None)
    # Prepare kwargs for 'bidirectional_best_hit' by replacing underscores with hyphens
    diamond_kwargs = {k.replace('_', '-'): v for k, v in kwargs.items()}

    # Paths for output files
    full_results_pkl_path = output / 'full_results.pkl'
    results_csv_path = output / 'AAI_table.csv'

    # Check if output files already exist
    if full_results_pkl_path.exists() and results_csv_path.exists():
        with open(full_results_pkl_path, 'rb') as file:
            full_results = pickle.load(file)
        result_df = pd.read_csv(results_csv_path)
        logging.warning(f"Output files '{full_results_pkl_path}' and '{results_csv_path}' already exist. Skipping computation.")
    else:
        # Compute the number of genes in the reference genome
        fasta_header_count = count_fasta_headers(reference)
        results = []
        full_results = {}

        # Perform BBH analysis in parallel
        with ProcessPoolExecutor(max_workers=threads) as executor:
            future_to_query = {
                executor.submit(
                    diamond_instance.bidirectional_best_hit,
                    fasta_file1=str(query),
                    fasta_file2=str(reference),
                    **diamond_kwargs
                ): query for query in queries
            }
            for future in tqdm(as_completed(future_to_query), total=len(queries), desc="Processing BBH", unit="query"):
                query = future_to_query[future]
                try:
                    bbh_df = future.result()
                    genome_id = Path(query).stem
                    homolog_count = len(bbh_df)
                    full_results[genome_id] = bbh_df
                    aai = bbh_df['pident'].mean()
                    results.append((genome_id, fasta_header_count, homolog_count, aai))
                except Exception as e:
                    logging.error(f"Query {query} generated an exception: {e}, {traceback.format_exc()}")
                    raise e

        # Save full_results to pickle file
        with open(full_results_pkl_path, 'wb') as file:
            pickle.dump(full_results, file)

        # Create result DataFrame and save to CSV
        result_df = pd.DataFrame(results, columns=['Genome_ID', 'Total_Genes_Count', 'Homologous_Genes_Count', 'AAI'])
        result_df['Percentage_Bidirectional_Best_Hits'] = (result_df['Homologous_Genes_Count'] / result_df['Total_Genes_Count']) * 100
        result_df.to_csv(results_csv_path, index=False)

    # Concatenate all BBH DataFrames and add 'genome_id' column
    all_bbh_data = pd.concat([df.assign(genome_id=genome_id) for genome_id, df in full_results.items()], ignore_index=True)

    gene_presence_count = all_bbh_data['sseqid'].value_counts()

    # Identify core genes based on core_threshold
    num_genomes = len(full_results)
    min_presence = max(int(num_genomes * core_threshold / 100), 1)  # Ensure at least one genome
    core_genes = gene_presence_count[gene_presence_count >= min_presence].index

    # Filter to core genes
    core_hits_data = all_bbh_data[all_bbh_data['sseqid'].isin(core_genes)]

    # Create core genome matrix
    coregenome_matrix = core_hits_data.pivot_table(
        index='sseqid', columns='genome_id', values='qseqid', aggfunc='first'
    ).reset_index()
    coregenome_matrix.rename(columns={'sseqid': 'Core_Gene'}, inplace=True)
    coregenome_matrix.columns.name = None  # Remove the columns name

    return result_df, coregenome_matrix

async def check_and_extract_gene_sequences(core_gene, coregenome_matrix, orf_fnas, output_dir):
    """
    Checks if alignment file exists, and if not, extracts sequences for a core gene from all genomes asynchronously.
    """
    output_file = Path(output_dir) / f"{core_gene}.aln"
    
    # Check if the alignment file already exists
    if output_file.exists():
        return core_gene, None, f"Alignment already exists for {core_gene}, skipping extraction..."

    sequences = []
    for genome_id in coregenome_matrix.columns[1:]:  # Excludes 'Core_Gene'
        gene_id = coregenome_matrix.loc[coregenome_matrix['Core_Gene'] == core_gene, genome_id].iloc[0]
        
        if pd.isna(gene_id):
            sequences.append(f">{genome_id}\n-\n")
        else:
            try:
                if genome_id in orf_fnas and gene_id in orf_fnas[genome_id]:
                    seq = orf_fnas[genome_id][gene_id]
                    sequences.append(f">{genome_id}\n{seq.seq}\n")
                else:
                    sequences.append(f">{genome_id}\n-\n")
            except Exception as e:
                logging.error(f"Could not process {genome_id} for {gene_id}: {e}")
                sequences.append(f">{genome_id}\n-\n")
        
        # Add a small delay to allow other tasks to run
        await asyncio.sleep(0)
    
    return core_gene, "".join(sequences), "Sequences extracted successfully"

def align_and_save_gene(core_gene, sequences, output_dir, mafft):
    """
    Aligns sequences for a core gene using MAFFT and saves the result.
    """
    output_file = Path(output_dir) / f"{core_gene}.aln"
    
    if sequences is None:  # This means the alignment file already exists
        return core_gene, f"Alignment already exists for {core_gene}, skipping..."
    
    if sequences.strip():  # Check if there are sequences to align
        start_time = time.time()
        mafft.align_and_save(sequences, str(output_file), algorithm="auto")
        align_time = time.time() - start_time
        return core_gene, f"Alignment saved for {core_gene} - Align time: {align_time:.2f}s" 
    else:
        return core_gene, f"No sequences to align for {core_gene}, skipping..."

async def extract_align_and_save_core_genes(coregenome_matrix, orf_fnas, output_dir, max_workers=4):
    """
    Asynchronously extracts core genome genes, aligns them concurrently as they become available,
    and frees memory for completed genes.

    Args:
        coregenome_matrix (pd.DataFrame): Matrix of core genes across genomes.
        orf_fnas (dict): Dictionary of ORF sequences for each genome.
        output_dir (str or Path): Directory to save aligned gene sequences.
        max_workers (int): Maximum number of concurrent workers for extraction and alignment.

    Returns:
        None
    """
    core_genes = coregenome_matrix['Core_Gene'].tolist()
    mafft = MAFFT()
    
    # Ensure output directory exists
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Create queues for extraction and alignment tasks
    extraction_queue = asyncio.Queue()
    alignment_queue = asyncio.Queue()

    # Count existing alignment files
    existing_alignments = sum(1 for gene in core_genes if (Path(output_dir) / f"{gene}.aln").exists())

    # Add extraction tasks to the queue
    for gene in core_genes:
        await extraction_queue.put(gene)
    
    # Add sentinel values to signal workers to terminate
    for _ in range(max_workers):
        await extraction_queue.put(None)
    
    # Initialize progress bars
    extraction_pbar = tqdm(total=len(core_genes), desc="Extracting genes", position=0)
    alignment_pbar = tqdm(total=len(core_genes), desc="Aligning genes", position=1)
    
    # Update extraction progress bar for existing alignments
    extraction_pbar.update(existing_alignments)
    alignment_pbar.update(existing_alignments)

    extraction_lock = asyncio.Lock()
    alignment_lock = asyncio.Lock()
    
    # Create a ThreadPoolExecutor for alignment tasks
    alignment_executor = ThreadPoolExecutor(max_workers=max_workers)
    
    async def extract_worker():
        while True:
            gene = await extraction_queue.get()
            if gene is None:
                # Signal alignment workers to terminate
                await alignment_queue.put(None)
                extraction_queue.task_done()
                break
            gene, sequences, message = await check_and_extract_gene_sequences(
                gene, coregenome_matrix, orf_fnas, output_dir
            )
            if sequences is not None:
                await alignment_queue.put((gene, sequences))
            else:
                logging.info(message)
            extraction_queue.task_done()
            async with extraction_lock:
                extraction_pbar.update(1)
    
    async def align_worker():
        while True:
            item = await alignment_queue.get()
            if item is None:
                alignment_queue.task_done()
                break
            gene, sequences = item
            try:
                result = await asyncio.get_event_loop().run_in_executor(
                    alignment_executor, align_and_save_gene, gene, sequences, output_dir, mafft
                )
                logging.info(result[1])
            except Exception as e:
                logging.error(f"Error processing gene {gene}: {e}")
            finally:
                alignment_queue.task_done()
                async with alignment_lock:
                    alignment_pbar.update(1)
    
    # Start extraction workers
    extraction_tasks = [asyncio.create_task(extract_worker()) for _ in range(max_workers)]
    
    # Start alignment workers
    alignment_tasks = [asyncio.create_task(align_worker()) for _ in range(max_workers)]
    
    # Wait for all extractions and alignments to complete
    await extraction_queue.join()
    await alignment_queue.join()
    
    # Close progress bars
    extraction_pbar.close()
    alignment_pbar.close()
    
    # Cancel worker tasks
    for task in extraction_tasks + alignment_tasks:
        task.cancel()
    
    # Shut down the executor
    alignment_executor.shutdown()

def compute_entropy(gene_file, core_genome_genes_dir):
    gene_id = os.path.splitext(gene_file)[0]
    aln_path = os.path.join(core_genome_genes_dir, gene_file)
    alignment = AlignIO.read(aln_path, 'fasta')
    entropy = calculate_shannon_entropy(alignment)
    return gene_id, entropy

def process_chunk(chunk_index, chunk_gene_files, entropy_dict, core_genome_genes_dir, entropy_threshold):
    try:
        # Temporary files for this chunk
        temp_file_all = tempfile.NamedTemporaryFile(mode='w+', delete=False, prefix=f'chunk_{chunk_index}_all_')
        temp_file_filtered = tempfile.NamedTemporaryFile(mode='w+', delete=False, prefix=f'chunk_{chunk_index}_filtered_')
        temp_files = {'all': temp_file_all, 'filtered': temp_file_filtered}

        concatenated_sequences_all = {}
        concatenated_sequences_filtered = {}
        genome_ids = None

        for gene_file in chunk_gene_files:
            gene_id = os.path.splitext(gene_file)[0]
            entropy = entropy_dict[gene_id]
            aln_path = os.path.join(core_genome_genes_dir, gene_file)
            alignment = AlignIO.read(aln_path, 'fasta')

            if genome_ids is None:
                genome_ids = [record.id for record in alignment]
                for genome_id in genome_ids:
                    concatenated_sequences_all[genome_id] = ''
                    concatenated_sequences_filtered[genome_id] = ''

            for record in alignment:
                genome_id = record.id
                seq = str(record.seq)
                concatenated_sequences_all[genome_id] += seq
                if entropy < entropy_threshold:
                    concatenated_sequences_filtered[genome_id] += seq

        # Write concatenated sequences to temporary files
        for genome_id in genome_ids:
            temp_files['all'].write(f">{genome_id}\n{concatenated_sequences_all[genome_id]}\n")
            temp_files['filtered'].write(f">{genome_id}\n{concatenated_sequences_filtered[genome_id]}\n")

        temp_files['all'].close()
        temp_files['filtered'].close()
        return temp_files['all'].name, temp_files['filtered'].name
    except Exception as e:
        logging.error(f"Error in chunk {chunk_index}: {e}")
        raise e

def plot_entropy_distribution(entropy_dict, threshold, output_dir):
    gene_ids = list(entropy_dict.keys())
    entropy_values = list(entropy_dict.values())
    
    plt.figure(figsize=(12, 6))
    
    # Scatter plot
    plt.scatter(range(len(gene_ids)), entropy_values, alpha=0.6, s=10)
    
    # Threshold line
    plt.axhline(y=threshold, color='r', linestyle='--', label=f'Threshold: {threshold:.4f}')
    
    # Interquartile range
    q1, q3 = np.percentile(entropy_values, [25, 75])
    iqr = q3 - q1
    plt.axhspan(q1, q3, alpha=0.2, color='green', label='Interquartile Range')
    
    # Median line
    median = np.median(entropy_values)
    plt.axhline(y=median, color='green', linestyle='-', label=f'Median: {median:.4f}')
    
    plt.xlabel('Genes')
    plt.ylabel('Entropy')
    plt.title('Distribution of Gene Entropies')
    plt.legend()
    
    # Adjust x-axis
    plt.xticks([])  # Remove x-axis ticks for clarity
    plt.xlim(0, len(gene_ids))
    
    # Add text for eliminated genes
    eliminated = sum(1 for e in entropy_values if e >= threshold)
    total = len(entropy_values)
    plt.text(0.02, 0.98, f"Eliminated: {eliminated}/{total} ({eliminated/total:.2%})",
             transform=plt.gca().transAxes, verticalalignment='top')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'entropy_distribution.png', dpi=300)
    plt.close()

def concatenate_alignments(core_genome_genes_dir, final_output, final_output_F, output_dir, max_workers=4):
    """
    Concatenate gene alignments into concatenated alignments for all genes and filtered genes in parallel.

    Parameters:
    - core_genome_genes_dir: Path to directory containing gene alignments
    - final_output: Path to output file for all genes concatenated alignment
    - final_output_F: Path to output file for filtered genes concatenated alignment
    - max_workers: Number of parallel workers to use
    """

    # Step 1: Collect entropy values and gene files
    entropy_dict = {}  # gene_id -> entropy_value
    gene_files = [f for f in os.listdir(core_genome_genes_dir) if f.endswith('.aln')]
    gene_files = sorted(gene_files)  # Ensure consistent order

    # Compute entropies in parallel
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(compute_entropy, gene_file, core_genome_genes_dir) for gene_file in gene_files]
        for future in tqdm(as_completed(futures), total=len(futures), desc="Computing entropies"):
            try:
                gene_id, entropy = future.result()
                entropy_dict[gene_id] = entropy
            except Exception as e:
                logging.error(f"Error computing entropy: {e}")
                raise e

    # Save entropy values to CSV
    entropy_df = pd.DataFrame.from_dict(entropy_dict, orient='index', columns=['Entropy'])
    entropy_df.index.name = 'Gene_ID'
    entropy_csv_path = output_dir / 'gene_entropies.csv'
    entropy_df.to_csv(entropy_csv_path)
    logging.info(f"Saved gene entropies to {entropy_csv_path}")

    # Step 2: Compute entropy threshold
    all_entropy = list(entropy_dict.values())
    q1, q3 = np.percentile(all_entropy, [25, 75])
    iqr = q3 - q1
    entropy_threshold = q3 + 1.5 * iqr

    # Count genes eliminated due to entropy
    genes_eliminated = sum(1 for entropy in all_entropy if entropy >= entropy_threshold)
    total_genes = len(all_entropy)
    
    # Log warning about eliminated genes
    logging.warning(f"{genes_eliminated} out of {total_genes} genes ({genes_eliminated/total_genes:.2%}) were eliminated due to high entropy (threshold: {entropy_threshold:.4f})")

    # Generate and save the entropy distribution plot
    plot_entropy_distribution(entropy_dict, entropy_threshold, output_dir)

    # Step 3: Divide genes into chunks
    chunk_size = max(1, len(gene_files) // max_workers)
    gene_chunks = [gene_files[i:i + chunk_size] for i in range(0, len(gene_files), chunk_size)]

    temp_file_names = []

    # Step 4: Process chunks in parallel
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for i, chunk in enumerate(gene_chunks):
            futures.append(executor.submit(
                process_chunk, i, chunk, entropy_dict, core_genome_genes_dir, entropy_threshold))
        for future in tqdm(as_completed(futures), total=len(futures), desc="Processing chunks"):
            result = future.result()
            temp_file_names.append(result)

    # Step 5: Combine temporary files
    def merge_temp_files(output_file, temp_file_paths):
        # Collect genome IDs to ensure consistent order
        genome_sequences = {}
        for temp_file in temp_file_paths:
            with open(temp_file, 'r') as f:
                for record in SeqIO.parse(f, 'fasta'):
                    genome_id = record.id
                    seq = str(record.seq)
                    if genome_id in genome_sequences:
                        genome_sequences[genome_id] += seq
                    else:
                        genome_sequences[genome_id] = seq
            os.unlink(temp_file)  # Delete temp file after reading

        # Write merged sequences to the output file
        with open(output_file, 'w') as f_out:
            for genome_id, seq in genome_sequences.items():
                f_out.write(f">{genome_id}\n{seq}\n")

    # Separate temp files for all genes and filtered genes
    temp_files_all = [pair[0] for pair in temp_file_names]
    temp_files_filtered = [pair[1] for pair in temp_file_names]

    # Merge temporary files into final outputs
    merge_temp_files(final_output, temp_files_all)
    merge_temp_files(final_output_F, temp_files_filtered)

    logging.info("Concatenation completed.")

def main():
    args, output_dir = parse_arguments()
    setup_logging(args.verbose, output_dir)
    
    input_files = [Path(x) for x in args.input]
    
    logging.info(f"Starting FABBIT version {__version__}")
    logging.debug(f"Input files: {', '.join(str(f) for f in input_files)}")
    logging.debug(f"Output directory: {output_dir}")

    subdirs = create_output_directories(output_dir)
    final_output = output_dir / 'core_genome.aln'
    final_output_F = output_dir / 'core_genome_F.aln'
    results_table_path = output_dir / 'AAI_table.csv'
    pyrodigal_orfs = output_dir / 'pyrodigal_orfs'
    core_genome_genes = output_dir / 'core_genome_genes'

    # Process genomes and predict ORFs
    process_genomes(input_files, pyrodigal_orfs, max_workers=args.threads)

    # Reading the orfs from the pyrodigal output
    orf_fnas = {}
    for fna in tqdm(list(pyrodigal_orfs.glob('*.fna')), desc="Reading ORFs"):
        orf_fnas[fna.stem] = pyfastx.Fasta(str(fna))

    diamond = DIAMOND()
    faa_sample_files = [pyrodigal_orfs / (x.stem + '.faa') for x in input_files]
    
    # Select reference genome
    if args.reference:
        reference_file = Path(args.reference)
        if reference_file not in input_files:
            raise FileNotFoundError(f"Reference file {reference_file} is not in the input files.")
        reference_faa = pyrodigal_orfs / (reference_file.stem + '.faa')
    else:
        reference_faa = faa_sample_files[0]
        logging.info(f"No reference genome specified. Using {reference_faa} as reference.")
        
    logging.info(f"Reference genome: {reference_faa}")
    logging.debug(f"Sample genomes: {', '.join(str(f) for f in faa_sample_files)}")
    
    results_table, coregenome_matrix = parallel_bidirectional_best_hit(
        reference=reference_faa,
        queries=faa_sample_files,
        diamond_instance=diamond,
        output=output_dir,
        core_threshold=args.core_threshold,
        threads=args.threads,
        sensitivity=args.sensitivity,
        evalue=args.evalue,
        query_cover=args.query_cover,
        max_target_seqs=args.max_target_seqs,
        id=args.id,
    )
    logging.info("Bidirectional best hit analysis completed")
    logging.debug(f"Results table:\n{results_table}")
    results_table.to_csv(results_table_path, index=False)
    coregenome_matrix.to_csv(output_dir / 'core_genes_matrix.csv', index=False)

    # Extract, align and save core genome genes
    asyncio.run(extract_align_and_save_core_genes(
        coregenome_matrix,
        orf_fnas,
        core_genome_genes,
        max_workers=args.threads
    ))

    # Now, concatenate the alignments
    concatenate_alignments(
        core_genome_genes_dir=core_genome_genes,
        final_output=final_output,
        final_output_F=final_output_F,
        output_dir=output_dir,
        max_workers=args.threads
    )

    logging.info("FABBIT processing completed")

if __name__ == "__main__":
    main()