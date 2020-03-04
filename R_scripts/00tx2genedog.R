#tx2gene for human, dog, mouse

library(tximport)
library(readr)
library(dplyr)
library(tidyverse)
setwd ("C:/Users/Alicia/Desktop/XSP2")

###Dogs
#find and download tabular file for canine from NCBI (feature_table.txt.gz)
download.file(url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/285/GCF_000002285.3_CanFam3.1/GCF_000002285.3_CanFam3.1_feature_table.txt.gz",
              destfile = "inputs/GCF_000002285.3_CanFam3.1_feature_table.txt.gz")
feat_table_cf <- read_tsv('inputs/GCF_000002285.3_CanFam3.1_feature_table.txt.gz')

#select relevant columns: product accession and symbol (in that order)
feat_table_cf <- select(feat_table_cf, product_accession, symbol)
feat_table_cf <- feat_table_cf%>%
  filter(!is.na(product_accession))%>%
  unique()

###relevant columns used to include GeneID, had to remove this since tximport only tolerates 2 colums. Will need to recall feat_tables and add back in for Ortho genes transformation

#make sure table is capturing all salmon counts, should be 82430 obs in canine
quant <- read_tsv("outputs/Canine_NK_Day_0_HVE5_2_15_19_quant.sf")
table(quant$Name %in% feat_table_cf$product_accession)

#output table to folder
write.table(feat_table_cf, "outputs/tx2gene/dogtx2gene.tsv", quote = F, row.names = F, sep = "\t")
