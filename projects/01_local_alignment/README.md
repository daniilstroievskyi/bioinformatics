# Project 1 - Local Alignment

## Introduction

Imagine that you are scientists working for the secret genetic laboratory GENE-X, which is engaged in research into mysterious biological samples found in a meteorite from space. 
After thorough analysis, you have discovered that the nucleotide and amino acid sequences in these samples are similar to those found in terrestrial organisms, but contain strange mutations.
Your task is to implement an algorithm that will allow you to **locally align these sequences** with existing proteins on Earth so that we can determine whether this is a known form of life or something completely new!
To accomplish this mission, you will use the **Smith-Waterman algorithm**, which uses **dynamic programming** to find the best local alignment of sequences. 
Get ready, because based on your results, it may be decided whether this is an extraordinary discovery or just contamination by terrestrial bacteria!


## Objectives

Implement a program that performs local alignment of two sequences using the Smith-Waterman algorithm. 
The program must be divided into two parts:

1. Alignment of nucleotide sequences: Implement local alignment of two nucleotide sequences from the input.

2. Alignment of amino acid sequences: Implement local alignment of two amino acid sequences using the _BLOSUM62_ or _PAM250_ substitution matrix.

The program must:

1. Read two sequences from input files in FASTA format.

2. Use dynamic programming to compute the optimal local alignment score and the corresponding alignment.

3. Correctly reconstruct the optimal local alignment from the computed score matrix.

4. Compare the results of your implementation with those obtained from existing tools.


## Submission and Grading Criteria

You will submit a single ZIP file containing:

- Source code of your implementation.

- Example input files (FASTA format) and output files showing the alignments.

- A brief PDF report (3 - 4 pages) describing your implementation, obtained results, and comparison with existing tools.


The project will be graded based on the following criteria:

1. DNA Sequence Alignment (**6 points**):

    - Load nucleotide sequences from FASTA files (**1 point**)

    - Correctly implement DP algorithm for nucleotide sequences (**2 points**)

    - Reconstruct optimal local alignment from score matrix (**2 points**)

    - Compare results with existing tools (**1 point**)

2. Amino Acid Sequence Alignment (**6 points**):

    - Load amino acid sequences from FASTA files (**1 point**)

    - Correctly implement DP algorithm for amino acid sequences using (**2 points**)

    - Use _BLOSUM62_ and _PAM250_ substitution matrix (**1 point**)

    - Reconstruct optimal local alignment from score matrix (**1 point**)

    - Compare results with existing tools (**1 point**)

3. Presentation and Documentation (**3 points**)

**Total: 15 points**

**Minimum passing score: 10 points (66%)**


> **Use of AI is strictly prohibited and will be treated as an academic misconduct!**
The project submission will be accompanied by AI-detection and short presentation of the work done in front of the instructor.
Failure to demonstrate that the work is your own will result in zero points for the project.


## Further Reading and Resources

- Smith TF, Waterman MS. Identification of common molecular subsequences. J Mol Biol. 1981 Mar 25;147(1):195-197. doi: [10.1016/0022-2836(81)90087-5](https://www.sciencedirect.com/science/article/abs/pii/0022283681900875). PMID: [7265238](https://pubmed.ncbi.nlm.nih.gov/7265238/).

- National Center for Biotechnology Information (NCBI) - BLAST: https://blast.ncbi.nlm.nih.gov/Blast.cgi

- EMBL-EBI - Sequence Alignment Tools: https://www.ebi.ac.uk/jdispatcher/psa

- EMBOSS WATER - Smith-Waterman local alignment tool: https://www.ebi.ac.uk/jdispatcher/psa/emboss_water
