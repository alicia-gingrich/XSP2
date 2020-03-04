#now join counts and feature tables for each species
#download tx2gene dog feature table again
download.file(url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/285/GCF_000002285.3_CanFam3.1/GCF_000002285.3_CanFam3.1_feature_table.txt.gz",
              destfile = "inputs/GCF_000002285.3_CanFam3.1_feature_table.txt.gz")
feat_table_cf <- read_tsv('inputs/GCF_000002285.3_CanFam3.1_feature_table.txt.gz')

#select relevant columns: product accession, symbol and GeneID (this time GeneID IS included!)
feat_table_cf <- select(feat_table_cf, product_accession, symbol, GeneID)
feat_table_cf <- feat_table_cf%>%
  filter(!is.na(product_accession))%>%
  unique()
head(feat_table_cf)

#check tables for overlap
table(all_ortho_entrez2$Dog %in% feat_table_cf$GeneID)
table(feat_table_cf$GeneID %in% all_ortho_entrez2$Dog)

#make GeneID column in feat_table_dog a character instead of numeric 
feat_table_cf$GeneID <- as.character(feat_table_cf$GeneID)

#join ortho_entrez and feat_table_dog to make a table with tx2gene and entrezID of 1:1 ortholog pairs
feat_table_cf_ortho <- left_join(feat_table_cf, all_ortho_entrez2, by = c("GeneID"="Dog"))
head(feat_table_cf_ortho)

#filter this new table to symbol, entrezID
feat_table_cf_ortho <- feat_table_cf_ortho%>%
  select(symbol, entrezID.x)%>%
  distinct()%>%
  filter(!is.na(entrezID.x))
head(feat_table_cf_ortho)
#should have 8137 obs. of 2 variables

#read in dog counts from local computer
dog_counts <- read.csv("C:/Users/Alicia/Desktop/XSP2/outputs/counts/dog_counts.csv")
head(dog_counts)

#join table with symbols and entrezID for all dog genes with human 1:1 orthos, to dog counts
dog_ortho_counts <- inner_join(dog_counts, feat_table_cf_ortho, by = c("X"="symbol"))
head(dog_ortho_counts)
#should have 8123 obs. of 10 variables

#save dog_ortho_counts to XSP2 (working directory) folder for DESeq later
write.csv(as.data.frame(dog_ortho_counts), 
          file="dog_ortho_counts.csv")
