###then humans
# read in file names 
files <- read_csv("samples_hs.csv")

# read in human counts
tx2gene <- read_tsv("outputs/tx2gene/humantx2gene.tsv")
human <- tximport(files = files$files, type = "salmon", tx2gene = tx2gene)
humancounts <- human$counts

# change rownames for provenance
colnames(humancounts) <- files$sample

# write to file
write.csv(humancounts, "outputs/counts/human_counts.csv", quote = F, row.names = T)
