library(readr)
library(dplyr)
download.file(url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_feature_table.txt.gz",
              destfile = "inputs/GCF_000001635.26_GRCm38.p6_feature_table.txt.gz")
feat_table_mm <- read_tsv('inputs/GCF_000001635.26_GRCm38.p6_feature_table.txt.gz')

#select relevant columns: product accession and symbol (in that order)
feat_table_mm <- select(feat_table_mm, product_accession, symbol)
feat_table_mm <- feat_table_mm%>%
  filter(!is.na(product_accession))%>%
  unique()

#make sure table is capturing all salmon counts, should be 120102 in mouse
quant <- read_tsv("outputs/quant_mm_nonribo/M1_R_quant/quant.sf")
table(quant$Name %in% feat_table_mm$product_accession)

#output table to folder
write.table(feat_table_mm, "outputs/tx2gene/mousetx2gene.tsv", quote = F, row.names = F, sep = "\t")
