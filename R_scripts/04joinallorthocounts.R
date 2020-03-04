###And finally, inner join all dog, human and mouse ortho counts tables
#read in human, dog, and mouse ORTHO counts from local computer if you need to
#dog_ortho_counts <- read.csv("C:/Users/Alicia/Desktop/XSP/outputs/counts/dog_ortho_counts.csv")
#human_ortho_counts <- read.csv("C:/Users/Alicia/Desktop/XSP/outputs/counts/human_ortho_counts.csv")
#mouse_ortho_counts <- read.csv("C:/Users/Alicia/Desktop/XSP/outputs/counts/mouse_ortho_counts.csv")

two_ortho_counts <- inner_join(dog_ortho_counts, human_ortho_counts, by = "entrezID.x")

all_ortho_counts <- inner_join(two_ortho_counts, mouse_ortho_counts, by = "entrezID.x")

#7013 obs or 26 variables

write.csv(as.data.frame(all_ortho_counts), 
          file="all_ortho_counts.csv")
