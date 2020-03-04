#find and download tabular file for human from NCBI(feature_table.txt.gz)
download.file(url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_feature_table.txt.gz",
              destfile = "inputs/GCF_000001405.39_GRCh38.p13_feature_table.txt.gz")
feat_table_hs <- read_tsv('inputs/GCF_000001405.39_GRCh38.p13_feature_table.txt.gz')

#select relevant columns: product accession and symbol (in that order)
feat_table_hs <- select(feat_table_hs, product_accession, symbol)
feat_table_hs <- feat_table_hs%>%
  filter(!is.na(product_accession))%>%
  unique()

#make sure table is capturing all salmon counts, should be 157862 in human
quant <- read_tsv("outputs/Human_flash_20180803_quant.sf")
table(quant$Name %in% feat_table_hs$product_accession)

#output table to folder
write.table(feat_table_hs, "outputs/tx2gene/humantx2gene.tsv", quote = F, row.names = F, sep = "\t")
