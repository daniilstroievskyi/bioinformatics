# Practicum 01: Linux Basics

## Introduction

In this practicum, we will cover the basics of using the Linux operating system. You will learn how to navigate the file system, manage files and directories, and use common command-line tools in the context of bioinformatics. The practicum is designed for beginners with little to no prior experience with Linux and environment management tools like Conda/Mamba.


## Objectives

By the end of this practicum, you will be able to:

1. Understand the Linux file system structure

2. Navigate the file system using command-line tools - `cd`, `ls`, `pwd`, `echo`, assign variables, etc.

3. Manage files and directories - `cp`, `mv`, `mkdir`, `rm`, `rmdir`, `touch`, `find`, etc.

4. Use basic command-line tools for file viewing - `cat`, `less`, `head`, `tail`, `grep`, etc.

5. Use command-line piping and redirection - `|`, `>`, `>>`, `<`

6. Understand file permissions and how to modify them - `chmod`, `chown`, `chgrp`

7. Compress and decompress files - `tar`, `gzip`, `gunzip`, etc.

8. Manipulate text files using command-line tools - `cut`, `sort`, `uniq`, `wc`, `awk`, `sed`, etc.

9. Use Linux-based text editors - `nano`, `vim`

10. Use Conda/Mamba for environment management and package installation


## Prerequisites

- A computer with Linux installed, access to a Linux server, or WSL (Windows Subsystem for Linux) set up on your Windows machine.

- Basic understanding of command-line interface concepts.

- Bioinfotmatics resources repository cloned to your local machine. If you haven't done this yet, follow the instructions in the main README file of the repository.


## Getting Started

If you are using Windows, ensure that you have WSL installed and set up. You can follow the official Microsoft guide [here](https://docs.microsoft.com/en-us/windows/wsl/install).
If you are using macOS or Linux, you can use the built-in terminal application.
Note that some commands may slightly differ between Linux distributions and macOS, but the core concepts remain the same.

Besides Linux, we will also use Conda or Mamba for environment management and package installation.
You don't have to install both, choose one of them, but Mamba is generally faster and more efficient.
It is also built on top of Conda, so if you install Mamba, you will also have Conda available.
You can follow the instructions on their respective websites:

- [Conda Installation Guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

- [Mamba Installation Guide](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html)

TL;DR you can install Mamba via Conda as follows:

```bash
# Download the Miniforge installer script
# macOS (Apple Silicon)
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Darwin-arm64.sh

# macOS (x86_64)
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Darwin-x86_64.sh

# Linux (x86_64)
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh

# Run the installer script
bash Miniforge3-*.sh

# Alternatively, you can download and install Miniconda from https://conda-forge.org/download/
```

All resources and files needed for this practicum are included in the bioinformatics resources repository. Navigate to the `bioinformatics/01_linux/` directory in your terminal to find the relevant materials.


## Further Reading and Resources

- [Linux Survival](https://linuxsurvival.com/) - An interactive tutorial to learn Linux commands

- [The Linux Command Line](http://linuxcommand.org/tlcl.php) - A comprehensive guide to the Linux command line

- [GNU Core Utilities](https://www.gnu.org/software/coreutils/manual/coreutils.html) - Official documentation for core Linux command-line tools

- [AWK Tutorial](https://www.grymoire.com/Unix/Awk.html) - Comprehensive guide to using AWK for text processing, useful for manipulating sequencing data files

- [Sed Tutorial](https://www.grymoire.com/Unix/Sed.html) - Comprehensive guide to using Sed for stream editing, useful for manipulating sequencing data files

- [Grep Tutorial](https://www.gnu.org/software/grep/manual/grep.html) - Official documentation for Grep, a powerful tool for searching text in files

- [Conda Documentation](https://docs.conda.io/projects/conda/en/latest/) - Official documentation for Conda

- [Mamba Documentation](https://mamba.readthedocs.io/en/latest/) - Official documentation for Mamba

- [Anaconda Repository](https://anaconda.org/) - A repository of Conda packages
