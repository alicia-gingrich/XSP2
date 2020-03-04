library(DESeq2)
setwd ("C:/Users/Alicia/Desktop/XSP2")

counts <- read.csv("outputs/counts/all_ortho_counts_entrezID.csv")
counts_unique <- distinct(counts, entrezID.x, .keep_all = TRUE)

write.csv(as.data.frame(counts_unique), 
          file="all_ortho_counts_entrezID.csv")

counts <- read.csv("outputs/counts/all_ortho_counts_entrezID.csv",row.names = 1)
head(counts)

#now try to apply the counts
counts <- apply(counts, 1:2, round)
samples <- read.csv("samples.csv")

dds<- DESeqDataSetFromMatrix(counts,
                             colData = samples,
                             design = ~ treatment+species + treatment:species)
dds <- DESeq(dds)

res1 <- results(dds, contrast=c("treatment","resting","coculture"))

res11 <- results(dds, contrast=c("species","human","canine"))
res12 <- results(dds, contrast=c("species","canine","mouse"))
res12 <- results(dds, contrast=c("species","mouse","human"))

resOrdered1 <- res1[order(res1$pvalue),] 
write.csv(as.data.frame(resOrdered1), 
          file="res1_treated_results.csv")

vsd <- vst(dds, blind=FALSE)

library(RColorBrewer)
pal <- colorRampPalette(brewer.pal(11, "RdYlGn"))(100)

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:100]
df <- as.data.frame(colData(dds)[,c("species","treatment")])
heatmap <- pheatmap(assay(vsd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
                    cluster_cols=TRUE, annotation_col=df, color = pal, legend = TRUE, annotation_colors = list(
                      species = c(human = "#FF00FF", canine = "#0000FF", mouse = "#660099"),
                      treatment = c(resting = "#CCFFFF", coculture = "#003399")
                    ))

plotPCA(vsd, intgroup=c("species", "treatment"))

library(ggplot2)

#customize PCA plot with ggplot function
pcaData <- plotPCA(vsd, intgroup=c("species", "treatment"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=species, shape=treatment)) +
  geom_point(size=3) +
  scale_color_manual(values = c(human = "#FF00FF", canine = "#0000FF", mouse = "#660099"))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  theme_light()




