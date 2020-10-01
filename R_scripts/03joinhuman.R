library(readr)
library(dplyr)
###here is the code for repeating in humans
#download tx2gene human feature table again
download.file(url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.38_GRCh38.p12/GCF_000001405.38_GRCh38.p12_feature_table.txt.gz",
              destfile = "inputs/GCF_000001405.39_GRCh38.p12_feature_table.txt.gz")
feat_table_hs <- read_tsv('inputs/GCF_000001405.39_GRCh38.p12_feature_table.txt.gz')

#select relevant columns: product accession, symbol and GeneID (this time GeneID IS included!)
feat_table_hs <- select(feat_table_hs, product_accession, symbol, GeneID)
feat_table_hs <- feat_table_hs%>%
  filter(!is.na(product_accession))%>%
  unique()
head(feat_table_hs)

#check tables for overlap
table(all_ortho_entrez2$Human %in% feat_table_hs$GeneID)
table(feat_table_hs$GeneID %in% all_ortho_entrez2$Human)

#make GeneID column in feat_table_hs a character instead of numeric 
feat_table_hs$GeneID <- as.character(feat_table_hs$GeneID)

#join all_ortho_entrez2 and feat_table_hs to make a table with tx2gene and entrezID of 1:1 ortholog pairs
feat_table_hs_ortho <- left_join(feat_table_hs, all_ortho_entrez2, by = c("GeneID"="Human"))
head(feat_table_hs_ortho)

#filter this new table to symbol, entrezID
feat_table_hs_ortho <- feat_table_hs_ortho%>%
  select(symbol, entrezID.x)%>%
  distinct()%>%
  filter(!is.na(entrezID.x))
head(feat_table_hs_ortho)
dim(feat_table_hs_ortho) #should have 12816 obs. of 2 variables

#read in human counts from local computer
human_counts <- read.csv("outputs/counts/human_counts.csv")
head(human_counts)

#join table with symbols and entrezID for all human genes with dog 1:1 orthos, to human counts
human_ortho_counts <- inner_join(human_counts, feat_table_hs_ortho, by = c("X"="symbol"))
head(human_ortho_counts) 
dim(human_ortho_counts) # should have 12804 obs. of 10 variables or 12614 of 10 vars

#save human_ortho_counts to XSP (working directory) folder for DESeq later
write.csv(as.data.frame(human_ortho_counts), 
          file="outputs/ortho_counts/human_ortho_counts.csv")
