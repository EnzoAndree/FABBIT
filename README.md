# 🧬 FABBIT: FAst coregenome alignment Based on Bidirectional best hIT 🚀

## 🌟 Overview

FABBIT is a powerful bioinformatics tool designed for fast core genome alignment based on bidirectional best hit analysis. It's particularly useful for comparative genomics studies, allowing researchers to identify and analyze core genes across multiple genomes efficiently.

## 🎉 Features

- 🧬 Predicts Open Reading Frames (ORFs) using Pyrodigal
- 💎 Performs Bidirectional Best Hit (BBH) analysis using DIAMOND
- 🔍 Identifies core genes across multiple genomes
- 🧩 Aligns core genes using MAFFT
- 🧮 Calculates Average Amino Acid Identity (AAI)
- 🔗 Generates concatenated core genome alignments
- 🎭 Filters genes based on Shannon entropy
- ⚡ Provides parallel processing capabilities for improved performance

## 🛠️ Installation

You can install Fabbit directly from GitHub using pip:

```bash
pip install git+https://github.com/EnzoAndree/FABBIT.git
```

## 🚀 Usage

After installation, you can use Fabbit from the command line:

```bash
fabbit -i input1.fasta input2.fasta -o output_directory -t 4
```

### 🎛️ Arguments

- `-i`, `--input`: Path to the input FASTA file(s) (required, multiple files allowed)
- `-o`, `--output`: Path to the output directory (required)
- `-t`, `--threads`: Number of threads to use (default: 1)
- `-r`, `--reference`: Path to reference genome file (optional)
- `--sensitivity`: DIAMOND sensitivity mode (default: "fast")
- `--evalue`: Maximum e-value to report alignments (default: 1e-6)
- `--query-cover`: Minimum query cover percentage (default: 95)
- `--max-target-seqs`: Maximum number of target sequences per query (default: 25)
- `--id`: Minimum identity percentage (default: 30)
- `--core-threshold`: Threshold for defining core genes (default: 95)
- `-v`, `--verbose`: Verbose level: 1=ERROR, 2=WARNING, 3=INFO, 4=DEBUG (default: 2)
- `-V`, `--version`: Show program's version number and exit

## 📊 Output

FABBIT generates the following main output files:

1. `core_genome_all.aln`: Concatenated alignment of all core genes
2. `core_genome_filtered.aln`: Concatenated alignment of core genes after entropy filtering
3. `AAI_table.csv`: Summary of Average Amino acid Identity (AAI) results
4. `core_genes_matrix.csv`: Matrix of core genes across genomes
5. `entropy_distribution.png`: Plot of gene entropy distribution
6. `gene_entropies.csv`: CSV file containing entropy values for each gene
7. `partition_all.txt`: Partition file for all core genes
8. `partition_filtered.txt`: Partition file for filtered core genes

## 📁 Directory Structure

The script will create the following subdirectories in the specified output directory:
- `pyrodigal_orfs`: Contains predicted ORFs
- `core_genome_genes`: Contains individual core gene alignments
- `logs`: Contains log files

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🐛 Troubleshooting

If you encounter any issues, please check the log file in the `logs` directory for detailed error messages. For further assistance, please open an issue on the GitHub repository.

## 🙏 Acknowledgements

FABBIT makes use of several open-source tools and libraries. We thank the developers of DIAMOND, MAFFT, Pyrodigal, and other dependencies for their valuable contributions to the scientific community.