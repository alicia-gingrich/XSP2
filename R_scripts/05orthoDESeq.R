library(tidyr)
library(dplyr)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)

# setwd ("C:/Users/Alicia/Desktop/XSP2")

counts <- read.csv("outputs/ortho_counts/all_ortho_counts.csv")
counts_unique <- counts %>%
  distinct(entrezID.x, .keep_all = TRUE) %>%
  transmute(entrezID = entrezID.x, across(everything())) %>%
  select(-X.1, -X.x, -X.y, -X, -entrezID.x)

write.csv(as.data.frame(counts_unique), 
          file="outputs/ortho_counts/all_ortho_counts_entrezID.csv", row.names = F, 
          quote = F)

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
                             design = ~ treatment+species + treatment:species)
dds <- DESeq(dds)

res1 <- results(dds, contrast=c("treatment", "nkp46_resting", "nkp46_activated"))

res11 <- results(dds, contrast=c("species", "human", "canine"))
dim(subset(res11, padj < .05))
res12 <- results(dds, contrast=c("species", "canine", "mouse"))
res12 <- results(dds, contrast=c("species", "mouse", "human"))
dim(subset(res12, padj < .05))



resOrdered1 <- res1[order(res1$pvalue),] 
write.csv(as.data.frame(resOrdered1), 
          file="outputs/deseq2/res1_treated_results.csv")

vsd <- vst(dds, blind=FALSE)


pal <- colorRampPalette(brewer.pal(11, "RdYlGn"))(100)


select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:100]
df <- as.data.frame(colData(dds)[ , c("species","treatment")])
heatmap <- pheatmap(assay(vsd)[select,], cluster_rows=TRUE, show_rownames=F,
                    cluster_cols=TRUE, annotation_col = df, color = pal, legend = TRUE, 
                    annotation_colors = list(
                      species = c(human = "#FF00FF", canine = "#0000FF", mouse = "#660099"),
                      treatment = c(nkp46_resting = "#CCFFFF", nkp46_activated = "#003399")
                    ))

plotPCA(vsd, intgroup=c("species", "treatment"))



#customize PCA plot with ggplot function
pcaData <- plotPCA(vsd, intgroup=c("species", "treatment"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=species, shape=treatment)) +
  geom_point(size=3) +
  scale_color_manual(values = c(human = "#FF00FF", canine = "#0000FF", mouse = "#660099"))+
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  coord_fixed() +
  theme_light()
