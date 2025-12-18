# Practicum 07: Annotation and Interpretation of Variants

## Introduction

This practicum will focus on the annotation and interpretation of genetic variants identified from sequencing data.
You will learn how to use tools such as ANNOVAR, GATK Funcotator, and VEP to annotate variants with functional information, population frequencies, and clinical significance.
Part of this practicum will also cover how to use public databases to interpret the potential impact of variants on health and disease.


## Objectives

By the end of this practicum, you will be able to:

1. Annotate variants using tools like GATK Funcotator and VEP.

2. Inspect the annotated variant files to extract relevant information.

3. Integrate functional annotations, population frequencies, and clinical significance into variant datasets.

4. Use public databases such as ClinVar, dbSNP, and gnomAD for variant interpretation.

5. Prioritize variants based on their potential impact on health and disease.

6. Generate comprehensive variant annotation reports for downstream analysis.


## Prerequisites

- Bioinformatics resources repository cloned to your local machine. If you haven't done this yet, follow the instructions in the main README file of the repository.

- VCF files containing variants called from cleaned BAM files from **Practicum 06: Variant Calling**.


## Getting Started

Additional software/tools required for this practicum:

- [Bcftools](http://www.htslib.org/doc/bcftools.html)

- [GATK](https://gatk.broadinstitute.org/hc/en-us)

- [maftools](https://bioconductor.org/packages/release/bioc/html/maftools.html)

GATK and Bcftools should be already installed, and we need to install maftools.
maftools is an R package, so make sure you have R and Bioconductor installed on your system.
You can install maftools by running the following commands in R:

```bash
wget https://cran.r-project.org/bin/linux/ubuntu/focal-cran40/README
sudo apt install r-base

R

# Inside R console, install BiocManager
install.packages("BiocManager")

# Install maftools using BiocManager
BiocManager::install("maftools")

# Verify installation
library(maftools)
```

In case you want to install R using conda, you can do it as follows:

```bash
conda create -n r_env -c bioconda r-base r-essentials bioconductor-maftools
conda activate r_env

R

# Inside R console, verify installation
library(maftools)
```

Alternatively, you can also use Docker to run R with maftools installed. Refer to [this guide](https://hub.docker.com/r/sevenbridges/maftools) for more details.

```bash
docker pull sevenbridges/maftools
docker run -it --name bioinf -v $(pwd):/data sevenbridges/maftools R

# Inside R console, verify installation
library(maftools)
```

Another thing you will need is the database files for variant annotation.
For GATK Funcotator, you can download the required database files by following the instructions [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035889931-Funcotator-Information-and-Tutorial).
You can download the pre-packaged database bundle using the GATK FuncotatorDataSourceDownloader tool:

```bash
gatk FuncotatorDataSourceDownloader \ 
   --somatic \ 
   --validate-integrity \ 
   --extract-after-download \ 
   --output path/to/funcotator_db_somatic
```

> Google has recently complicated access to the resource files hosted on Google Cloud, so the above command may not work as expected.
If this doesn't work, contact the course instructor for alternative download links.

## Further Reading and Resources

- [Samtools Documentation](http://www.htslib.org/doc/samtools.html)

- [Qualimap Documentation](http://qualimap.bioinfo.cipf.es/doc_html/)

- [GATK Documentation](https://gatk.broadinstitute.org/hc/en-us)



