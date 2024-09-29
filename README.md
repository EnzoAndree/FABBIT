# Fabbit

FABBIT: FAst coregenome alignment Based on Bidirectional best hIT

## Overview

FABBIT is a powerful bioinformatics tool designed for fast core genome alignment based on bidirectional best hit analysis. It's particularly useful for comparative genomics studies, allowing researchers to identify and analyze core genes across multiple genomes efficiently.

## Features

- Predicts Open Reading Frames (ORFs) using Pyrodigal
- Performs Bidirectional Best Hit (BBH) analysis using DIAMOND
- Identifies core genes across multiple genomes
- Aligns core genes using MAFFT
- Calculates Average Amino Acid Identity (AAI)
- Generates concatenated core genome alignments
- Filters genes based on Shannon entropy
- Provides parallel processing capabilities for improved performance

## Installation

You can install Fabbit directly from GitHub using pip:

```bash
pip install git+https://github.com/EnzoAndree/FABBIT.git
```

## Usage

After installation, you can use Fabbit from the command line:

```bash
fabbit -i input1.fasta input2.fasta -o output_directory -t 4
```

Arguments:
- `-i`, `--input`: Path to the input FASTA file(s) (required)
- `-o`, `--output`: Path to the output directory (required)
- `-t`, `--threads`: Number of threads to use (default: 1)

The script will create the following subdirectories in the specified output directory:
- `pyrodigal_orfs`
- `core_genome_genes`
- `core_genome_concatenated`

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.