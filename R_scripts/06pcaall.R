library(dplyr)
library(readr)
library(ggplot2)
library(DESeq2)
library(tibble)
counts <- read.csv("outputs/ortho_counts/all_ortho_counts_entrezID.csv", row.names = 1)

#now try to apply the counts
counts <- apply(counts, 1:2, round)
samples <- read.csv("samples_test_lab_swap.csv") 
samples$treatment <- factor(samples$treatment, levels =c("nkp46_activated", "nkp46_resting", "cd5_resting"))
counts <- counts[ , order(match(colnames(counts), samples$sample))]

vsd_all <- vst(counts)
pca_all <- prcomp(t(vsd_all))
# calculate percent variance explained
eigs <- pca_all$sdev^2
pca1 <- round((eigs[1] / sum(eigs))*100, digits = 0)
pca2 <-  round((eigs[2] / sum(eigs))*100, digits = 0)

pca_all <- pca_all$x %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  left_join(samples)

ggplot(pca_all, aes(PC1, PC2, color=species, shape=treatment)) +
  geom_point(size=3) +
  scale_color_manual(values = c(human = "#FF00FF", canine = "#0000FF", mouse = "#660099"))+
  xlab(paste0("PC1: ", pca1, "% variance")) +
  ylab(paste0("PC2: ", pca2, "% variance")) + 
  coord_fixed() +
  theme_light() 
