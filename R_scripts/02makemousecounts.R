###then mice
# read in file names 
files <- read_csv("samples_mm.csv")

# read in mouse counts
tx2gene <- read_tsv("outputs/tx2gene/mousetx2gene.tsv")
mouse <- tximport(files = files$files, type = "salmon", tx2gene = tx2gene)
mousecounts <- mouse$counts

# change rownames for provenance
colnames(mousecounts) <- files$sample

# write to file
write.csv(mousecounts, "outputs/counts/mouse_counts.csv", quote = F, row.names = T)
