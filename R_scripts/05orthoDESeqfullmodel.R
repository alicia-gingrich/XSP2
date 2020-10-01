library(tidyr)
library(dplyr)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)

# setwd ("C:/Users/Alicia/Desktop/XSP2")

# counts <- read.csv("outputs/ortho_counts/all_ortho_counts.csv")
# counts_unique <- counts %>%
#   distinct(entrezID.x, .keep_all = TRUE) %>%
#   transmute(entrezID = entrezID.x, across(everything())) %>%
#   select(-X.1, -X.x, -X.y, -X, -entrezID.x)
# 
# write.csv(as.data.frame(counts_unique), 
#           file="outputs/ortho_counts/all_ortho_counts_entrezID.csv", row.names = F, 
#           quote = F)

counts <- read.csv("outputs/ortho_counts/all_ortho_counts_entrezID.csv", row.names = 1)
head(counts)

#now try to apply the counts
counts <- apply(counts, 1:2, round)
samples <- read.csv("samples.csv") %>%
  filter(treatment != "cd5_resting")



counts <- counts[ , colnames(counts) %in% samples$sample]
counts <- counts[ , order(match(colnames(counts), samples$sample))]
dds<- DESeqDataSetFromMatrix(counts,
                             colData = samples,
                             design = ~ species_treatment)
dds <- DESeq(dds)

mouse_human_resting <- results(dds, contrast=c("species_treatment", "mouse_nkp46_resting", "human_nkp46_resting"))
dim(subset(mouse_human_resting, padj < .05))
write_tsv(mouse_human_resting%>%
            as.data.frame() %>%
            rownames_to_column("entrezID"), "outputs/deseq2/mouse_human_resting.tsv")
mouse_human_activated <- results(dds, contrast=c("species_treatment", "mouse_nkp46_activated", "human_nkp46_activated"))
dim(subset(mouse_human_activated, padj < .05))
write_tsv(mouse_human_activated %>%
            as.data.frame() %>%
            rownames_to_column("entrezID"), "outputs/deseq2/mouse_human_activated.tsv")

canine_human_resting <- results(dds, contrast=c("species_treatment", "canine_nkp46_resting", "human_nkp46_resting"))
dim(subset(canine_human_resting, padj < .05))
write_tsv(canine_human_resting %>%
            as.data.frame() %>%
            rownames_to_column("entrezID"), "outputs/deseq2/canine_human_resting.tsv")

canine_human_activated <- results(dds, contrast=c("species_treatment", "canine_nkp46_activated", "human_nkp46_activated"))
dim(subset(canine_human_activated, padj < .05))
write_tsv(canine_human_activated %>%
            as.data.frame() %>%
            rownames_to_column("entrezID"), "outputs/deseq2/canine_human_activated.tsv")



