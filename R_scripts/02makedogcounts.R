#make counts for human, dog and mouse on their own

###first dog
# read in file names 
files <- read_csv("samples_cf.csv")

# read in dog counts
tx2gene <- read_tsv("outputs/tx2gene/dogtx2gene.tsv")
canine <- tximport(files = files$files, type = "salmon", tx2gene = tx2gene)
caninecounts <- canine$counts

# change rownames for provenance
colnames(caninecounts) <- files$sample

# write to file
write.csv(caninecounts, "outputs/counts/dog_counts.csv", quote = F, row.names = T)
