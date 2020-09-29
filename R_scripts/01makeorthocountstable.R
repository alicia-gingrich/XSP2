#Make Ortho Counts for Human, Dog, Mouse

#read in Human and Dog Orthologous genes list from OMA with Entrez notation 
ortho_entrez1 <- read_tsv("Ortho_HSA_CFA_Entrez.txt",
                          col_names = c("Human", "Dog", "Pair", "entrezID"))
head(ortho_entrez1)

#filter ortho_entrez1 for 1:1 matches 
ortho_entrez1 <- filter(ortho_entrez1, Pair=="1:1")
table(is.na(ortho_entrez1$entrezID))
head(ortho_entrez1)

ortho_entrez3 <- read_tsv("Ortho_CFA_MUS_Entrez.txt",
                          col_names = c("Dog", "Mouse", "Pair", "entrezID"))
head(ortho_entrez3)

#filter ortho_entrez3 for 1:1 matches 
ortho_entrez3 <- filter(ortho_entrez3, Pair=="1:1")
table(is.na(ortho_entrez3$entrezID))
head(ortho_entrez3)

#time to join the human, dog and mouse orthologs 
all_ortho_entrez <- inner_join(ortho_entrez1, ortho_entrez3, by = c("Dog"="Dog"))
head(all_ortho_entrez)

#remove the second Pair.y, entrezID.y column
all_ortho_entrez2 <- select(all_ortho_entrez, -(Pair.y))
all_ortho_entrez2 <- select(all_ortho_entrez2, -(entrezID.y))
head (all_ortho_entrez2)
#14278 obs of 5 variables (Human, dog, Pair, entrezID, mouse)

write.table(all_ortho_entrez2, "outputs/all_ortho_entrez2.tsv", quote = F, row.names = F, sep = "\t")
