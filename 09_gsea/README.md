# Practicum 09: Gene Set Enrichment Analysis (GSEA)

## Introduction

In this practicum, we will explore Gene Set Enrichment Analysis (GSEA), a powerful method for interpreting gene expression data by determining whether predefined sets of genes show statistically significant differences between two biological states.
You will learn how to perform GSEA using R and Bioconductor packages, visualize the results, interpret the biological significance of enriched gene sets, and use web-based tools to investigate metabolic and signaling pathways associated with differentially expressed genes.


## Objectives

By the end of this practicum, you will be able to:

1. Understand the principles and applications of Gene Set Enrichment Analysis (GSEA).

2. Perform GSEA using R and Bioconductor packages.

3. Visualize and interpret GSEA results.

4. Utilize web-based tools to explore metabolic and signaling pathways related to differentially expressed genes.


## Prerequisites

- Bioinformatics resources repository cloned to your local machine. If you haven't done this yet, follow the instructions in the main README file of the repository.

- Results from the RNA-Seq practicum, including the ranked list of differentially expressed genes based on log2 fold change and p-values.


## Getting Started

Additional software/tools required for this practicum:

- [Salmon](https://salmon.readthedocs.io/en/latest/salmon.html)

- [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)

- [enrichplot](https://bioconductor.org/packages/release/bioc/html/enrichplot.html)



You should have Salmon and DESeq2 installed from the RNA-Seq practicum.
We only need to install R packages for GSEA analysis.

```bash
R

# Inside R console, install BiocManager
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "DOSE", "enrichplot", "biomaRt", "tidyverse", "ggtree", "pathview"))

# Verify installation
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(biomaRt)
library(tidyverse)
library(ggtree)
library(pathview)
```

In case you want to install R using conda, you can do it as follows:

```bash
conda install -c bioconda bioconductor-clusterprofiler bioconductor-org.hs.eg.db bioconductor-dose bioconductor-enrichplot bioconductor-biomart r-tidyverse bioconductor-pathview
```

R

# Inside R console, verify installation
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(biomaRt)
library(tidyverse)
# Note: ggtree is not available in bioconda, so you may need to install it separately using BiocManager
```


## Further Reading and Resources

- [Biomedical Knowledge Mining Book](https://yulab-smu.top/biomedical-knowledge-mining-book/)

- [KEGG Pathway Database](https://www.genome.jp/kegg/pathway.html)

- [Reactome Pathway Database](https://reactome.org/)

- [Panther Pathway Database](https://www.pantherdb.org/)

