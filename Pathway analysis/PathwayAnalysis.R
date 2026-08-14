if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("tidyverse")
BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(pathview)
files=c("edgeR_input.csv")
Pseudo=read.csv("PseudomonasPAO1.tsv", sep="\t", header=TRUE)
map=Pseudo[,c("Symbol","Gene.ID")]



for(i in 1:length(files)){
  # Formating
  input=read.csv(files[i], header=TRUE)
  df=input
  dubs=df$gene_name[duplicated(df$gene_name)]
  df=df[!(df$gene_name %in% dubs),]
  
  symbols_df <- data.frame(Symbol = df$gene_name, logFC=df$logFC)
  matched_data <- merge(symbols_df, map, by = "Symbol", all.x = TRUE)
  
  original_gene_list = setNames(matched_data$logFC, matched_data$gene_name) %>% na.omit() %>% sort(decreasing = TRUE)
  kegg_gene_list = setNames(matched_data$logFC, matched_data$Gene.ID) %>% na.omit() %>% sort(decreasing = TRUE)
  
  # Test KEGG pathways (~84)
  
  kegg_organism = "pae"
  kk2 <- gseKEGG(geneList     = kegg_gene_list,
                 organism     = kegg_organism,
                 nPerm        = 10000,
                 minGSSize    = 3,
                 maxGSSize    = 800,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "none",
                 keyType       = "ncbi-geneid")
  write.csv(kk2, paste0("KEGG_new", files[i]))
  
  
  # Colors pathway pictures
  
  for(j in 1:nrow(kk2)){
    if(kk2[j, "Description"] != "Metabolic pathways") {
      
      try({
        pathview(
          gene.data = kegg_gene_list,
          pathway.id = kk2[j, "ID"],
          species = "pae",
          out.suffix = kk2[j, "Description"],
          out.dir = file_path)
      })
    }
  }
  
  
}

