###here is the code for repeating in mice
#download tx2gene mice feature table again
download.file(url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_feature_table.txt.gz",
              destfile = "inputs/GCF_000001635.26_GRCm38.p6_feature_table.txt.gz")
feat_table_mm <- read_tsv('inputs/GCF_000001635.26_GRCm38.p6_feature_table.txt.gz')

#select relevant columns: product accession, symbol and GeneID (this time GeneID IS included!)
feat_table_mm <- select(feat_table_mm, product_accession, symbol, GeneID)
feat_table_mm <- feat_table_mm%>%
  filter(!is.na(product_accession))%>%
  unique()
head(feat_table_mm)

#check tables for overlap
table(all_ortho_entrez2$Mouse %in% feat_table_mm$GeneID)
table(feat_table_mm$GeneID %in% all_ortho_entrez2$Mouse)

#make GeneID column in feat_table_hs a character instead of numeric 
feat_table_mm$GeneID <- as.character(feat_table_mm$GeneID)

#join all_ortho_entrez2 and feat_table_mm to make a table with tx2gene and entrezID of 1:1 ortholog pairs
feat_table_mm_ortho <- left_join(feat_table_mm, all_ortho_entrez2, by = c("GeneID"="Mouse"))
head(feat_table_mm_ortho)

#filter this new table to symbol, entrezID
feat_table_mm_ortho <- feat_table_mm_ortho%>%
  select(symbol, entrezID.x)%>%
  distinct()%>%
  filter(!is.na(entrezID.x))
head(feat_table_mm_ortho)
#should have 12473 obs. of 2 variables

#read in mouse counts from local computer
mouse_counts <- read.csv("C:/Users/Alicia/Desktop/XSP2/outputs/counts/mouse_counts.csv")
head(mouse_counts)

#join table with symbols and entrezID for all human genes with dog 1:1 orthos, to human counts
mouse_ortho_counts <- inner_join(mouse_counts, feat_table_mm_ortho, by = c("X"="symbol"))
head(mouse_ortho_counts)
#should have 12461 obs. of 10 variables

#save mouse_ortho_counts to XSP (working directory) folder for DESeq later
write.csv(as.data.frame(mouse_ortho_counts), 
          file="mouse_ortho_counts.csv")







