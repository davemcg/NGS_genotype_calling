
args <- commandArgs(trailingOnly=TRUE)
#When testing, comment out line above and use the line below.
#args <- c("Y:/crestAnno/forR.txt",
#          "Y:/resources/CRESTvariantClassification.xlsx", "Y:/crestAnno/annoted.xlsx")

library(tidyverse)
library(readxl)

annovar <- read_tsv(args[1], col_names = TRUE, na = c("NA", "", "None", "."), col_types = cols(.default = col_character())) %>%
  type_convert()

classificationDF <- read_xlsx(args[2], sheet = "SV", na = c("NA", "", "None", ".")) %>%
  select("leftChr", "leftPos", "leftStrand", "rightChr", "rightPos", "rightStrand", "SVtype",
         "classification", "popAF", "note")
classificationDF$classification = factor(classificationDF$classification, levels = c("Pathogenic", "Likely pathogenic", "VOUS", "Not classified", "Likely benign", "Benign", "Artifact")) 

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

annovar$leftAA <- sapply(1:nrow(annovar), function(x) {mapHGMD(annovar[x, "leftAA"])})
annovar$leftIntronic <- sapply(1:nrow(annovar), function(x) {mapHGMD(annovar[x, "leftIntronic"])})
annovar$rightAA <- sapply(1:nrow(annovar), function(x) {mapHGMD(annovar[x, "rightAA"])})
annovar$rightIntronic <- sapply(1:nrow(annovar), function(x) {mapHGMD(annovar[x, "rightIntronic"])})

#"leftGene", "leftIntronic", "leftAA", "rightGene", "rightIntronic", "rightAA",
annoted <- left_join(annovar, classificationDF, by = c("leftChr", "leftPos", "leftStrand", "rightChr", "rightPos", "rightStrand", "SVtype")) %>% 
  replace_na(list(classification = "Not classified")) %>% 
  mutate(consensusSequences = sub("^","https://genome.ucsc.edu/cgi-bin/hgBlat?type=BLAT%27s+guess&userSeq=", consensusSequences))

annoted1 <- left_join(annoted, panelGene, by = c("leftGene" = "gene")) %>% rename(leftclass = panel_class) %>% 
  replace_na(list(leftclass = "Other"))
annoted2 <- left_join(annoted1, panelGene, by = c("rightGene" = "gene")) %>%  rename(rightclass = panel_class) %>% 
  replace_na(list(rightclass = "Other")) %>% 
  mutate(eyeGene = ifelse(leftclass %in% c("Dx", "Candidate") | rightclass %in% c("Dx", "Candidate"), 1, 0)) %>% 
  select(sample, eyeGene, sumReads, indelSize, leftChr, leftPos, leftStrand, NumofLeftSoftClippedReads, 
         rightChr, rightPos, rightStrand, NumofRightSoftClippedReads, SVtype, 
         leftGene, leftIntronic, leftAA, rightGene, rightIntronic, rightAA, classification, popAF, note, consensusSequences, everything()) %>% 
  arrange(classification, desc(eyeGene), desc(sumReads))

openxlsx::write.xlsx(annoted2, file = args[5])