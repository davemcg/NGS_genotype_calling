
args <- commandArgs(trailingOnly=TRUE)
#When testing, comment out line above and use the line below.
# args <- c("W:/ddl_nisc_custom_capture/042020/scramble/107910Y.forR.txt",
#           "Y:/resources/SCRAMBLEvariantClassification.xlsx", "Z:/OGL_NGS/variant_prioritization/data/OGLv1_panel_DxORcandidate.tsv",
#           "W:/ddl_nisc_custom_capture/042020/scramble/107910Y.annoted.xlsx")

library(tidyverse)
library(readxl)

annovar <- read_tsv(args[1], col_names = TRUE, na = c("NA", "", "None", "."), col_types = cols(.default = col_character())) %>% 
  mutate(sample = args[6]) %>% 
  mutate(Intronic = gsub("0>-", "", Intronic))
# Do not use type_convert() because it converts chr:pos.

classificationDF <- read_xlsx(args[2], sheet = "Variant", na = c("NA", "", "None", ".")) %>% 
  select("Insertion", "MEI_Family", "Insertion_Direction", "classification", "popAF", "note")

#panelGene <- read_tsv(args[3], col_names = TRUE, col_types = cols(.default = col_character())) %>% select(gene, panel_class)
panelGene <- read_xlsx(args[3], sheet = "analysis", na = c("NA", "", "None", ".")) %>% select(gene, panel_class)

HGMD <- read_tsv(args[4], col_names = TRUE, na = c("NA", "", "None", "."), col_types = cols(.default = col_character()))
hgmdNM <- dplyr::pull(HGMD, name)

mapHGMD <- function(x){
  x = as.character(x)
  annovarAA <- as.list(strsplit(x, ",|;"))[[1]] 
  if (length(annovarAA) > 1) {
    hgmdVar <- purrr::keep(annovarAA, str_extract_all(annovarAA, "NM_[[:digit:]]+") %in% hgmdNM) 
    if (length(hgmdVar) == 0) {
      return(x)
    } else if (length(hgmdVar) == 1) { return(hgmdVar) }
    else { return(paste(hgmdVar, collapse = ",")) }
  } else {
    return(x)
  }
}

annovar$AA <- sapply(1:nrow(annovar), function(x) {mapHGMD(annovar[x, "AA"])})
annovar$Intronic <- sapply(1:nrow(annovar), function(x) {mapHGMD(annovar[x, "Intronic"])})

#annovar <- annovar %>% unite("refgenewithver", Intronic, AA, sep = ",", remove = TRUE, na.rm = TRUE)

annoted <- left_join(annovar, classificationDF, by = c("Insertion", "MEI_Family", "Insertion_Direction")) %>% 
  replace_na(list(classification = "Not classified")) 

annoted$classification = factor(annoted$classification, levels = c("Pathogenic", "Likely pathogenic", "VOUS", "Not classified", "Likely benign", "Benign", "Artifact")) 

annoted1 <- left_join(annoted, panelGene, by = c("Gene" = "gene")) %>%
  replace_na(list(panel_class = "Other")) %>% 
  mutate(eyeGene = ifelse(panel_class %in% c("Dx", "Candidate"), 1, 0)) %>% 
  mutate(Clipped_Reads_In_Cluster = as.integer(Clipped_Reads_In_Cluster)) %>% 
  select(eyeGene, Insertion, MEI_Family, Insertion_Direction, Clipped_Reads_In_Cluster, Alignment_Score, 
         Alignment_Percent_Length, Alignment_Percent_Identity, Clipped_Sequence, Clipped_Side, Start_In_MEI, Stop_In_MEI, 
         polyA_Position, polyA_Seq, polyA_SupportingReads, TSD, TSD_length, panel_class, Func_refGene, Gene, Intronic, AA, classification, popAF, sample, note) %>% 
  filter(classification %in% c("Pathogenic", "Likely pathogenic", "VOUS", "Not classified"), !grepl("GL", Insertion), Func_refGene %in% c("splicing", "exonic", "UTR5", "UTR3", "upstream")) %>% 
  arrange(classification, desc(eyeGene), desc(Clipped_Reads_In_Cluster))

if (dim(annoted1)[1] == 0) {
  write_tsv(annoted1, path = args[5])
} else {
  openxlsx::write.xlsx(annoted1, file = args[5])
}
