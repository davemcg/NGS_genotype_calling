
args <- commandArgs(trailingOnly=TRUE)
#When testing, comment out line above and use the line below.
# args <- c("W:/ddl_nisc_custom_capture/042020/scramble/107910Y.forR.txt",
#           "Y:/resources/SCRAMBLEvariantClassification.xlsx", "Z:/OGL_NGS/variant_prioritization/data/OGLv1_panel_DxORcandidate.tsv",
#           "W:/ddl_nisc_custom_capture/042020/scramble/107910Y.annoted.xlsx")

library(tidyverse)
library(readxl)

annoted <- read_tsv(args[1], col_names = TRUE, na = c("NA", "", "None", "."), col_types = cols(.default = col_character())) %>% 
  mutate(sample = args[6]) %>% 
  mutate(eyeGene = as.integer(eyeGene), Clipped_Reads_In_Cluster = as.integer(Clipped_Reads_In_Cluster))
#do not use type_convert() because it converts chr:pos.

annoted$classification = factor(annoted$classification, levels = c("Pathogenic", "Likely pathogenic", "VOUS", "Not classified", "Likely benign", "Benign", "Artifact")) 

annoted1 <- annoted %>%
  arrange(classification, desc(eyeGene), desc(Clipped_Reads_In_Cluster))

if (dim(annoted1)[1] == 0) {
  write_tsv(annoted1, path = args[2])
} else {
  openxlsx::write.xlsx(annoted1, file = args[2])
}