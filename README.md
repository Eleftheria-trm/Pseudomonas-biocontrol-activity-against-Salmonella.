
#Figures_Data_code

This folder contains all raw data for Figures 1b, 2, 4b and 4c. 
PK stands for PK307-4, and STM50 for Salmonella enterica Typhimurium 14028S tagged with lacIZ/kan. "Blue" corresponds to measured Salmonella colonies and "white" to measured Pseudomonas colonies. 

# Pseudomonas RNA seq-deseq folder

This repository contains R scripts and data for RNA sequencing Differential expression analysis on Pseudomonas in response to Salmonella presence in vitro and in planta
, supporting the research findings published in the companion paper (see citation)

## Directory Contents


### :file_folder: [deseq](./deseq/README.md)

Contains RNA-seq differential expression analysis.

## Usage

Run R scripts individually:

```bash
Rscript deseq/ReportForPaperFinal.R
```

## Requirements

R packages: `readxl`, `dplyr`, `tidyr`, `ggplot2`, `DESeq2`, `data.table`, `emmeans`, `lme4`, `ggtext`, `patchwork`, `openxlsx`

Note that a `./utils/install_dependencies.R` script can help install these.

## Citation

Vimont, N., Bastkowski, S., Savva, G. M., Bloomfield, S. J., Mather, A. E., Webber, M. A., & Trampari, E. (2025).  
*Environmental context as a key driver of Pseudomonas' biocontrol activity against Salmonella*.  
**bioRxiv**. [https://doi.org/10.1101/2025.06.23.661019](https://doi.org/10.1101/2025.06.23.661019)  
[Full text PDF](https://www.biorxiv.org/content/early/2025/06/23/2025.06.23.661019.full.pdf)  
[View at bioRxiv](https://www.biorxiv.org/content/early/2025/06/23/2025.06.23.661019)  


