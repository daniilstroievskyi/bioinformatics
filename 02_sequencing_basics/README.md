# Practicum 02: Basic Operations with Sequencing Data

## Introduction

In this practicum, we will focus on basic operations with sequencing data files commonly used in bioinformatics, such as FASTQ and FASTA files. 
You will learn how to manipulate these files using command-line tools, including viewing, filtering, and extracting relevant information. 
You will also practice working with tools that are essential for handling sequencing data in bioinformatics workflows.


## Objectives

By the end of this practicum, you will be able to:

1. Understand the structure and format of FASTQ and FASTA files.

2. View and inspect sequencing data files using command-line tools - `head`, `tail`, `less`, `cat`, etc.

3. Extract specific sequences or reads from FASTQ/FASTA files using tools like `grep`, `sed`, and `awk`.

4. Filter sequencing data based on quality scores or sequence content.

5. Convert between FASTQ and FASTA formats.

6. Use command-line tools to summarize sequencing data (e.g., counting reads, calculating average read length).


## Prerequisites

- Bioinformatics resources repository cloned to your local machine. If you haven't done this yet, follow the instructions in the main README file of the repository.


## Getting Started

Additional software/tools required for this practicum:

- [Seqtk](https://github.com/lh3/seqtk)

- [SeqKit](https://github.com/shenwei356/seqkit)

- [Samtools](http://www.htslib.org/)

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

You can install these tools using `conda/mamba`:

```bash
conda install -c bioconda seqtk seqkit samtools fastqc
```

All resources and files needed for this practicum are included in the bioinformatics resources repository. Navigate to the `bioinformatics/02_sequencing_basics/` directory in your terminal to find the relevant materials.


## Further Reading and Resources

- [FASTQ Format Specification](https://help.basespace.illumina.com/files-used-by-basespace/fastq-files) - Overview of the FASTQ file format

- [FASTA Format Specification](https://en.wikipedia.org/wiki/FASTA_format) - Overview of the FASTA file format

- [FastQC Manual](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/) - Documentation for FastQC tool with explanations of plots
