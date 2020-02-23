
args <- commandArgs(trailingOnly=TRUE)
#When testing, comment out line above and use the line below.
#args <- c("Y:/crestAnno/forR.txt",
#          "Y:/resources/CRESTvariantClassification.xlsx", "Y:/crestAnno/annoted.xlsx")

library(tidyverse)
library(readxl)

predSV <- read_tsv(args[1], col_names = TRUE, na = c("NA", "", "None", "."), col_types = cols(.default = col_character())) %>%
  type_convert()

classificationDF <- read_xlsx(args[2], sheet = "SV", na = c("NA", "", "None", ".")) %>%
  select("leftChr", "leftPos", "leftStrand", "rightChr", "rightPos", "rightStrand", "SVtype",
         "classification", "popAF", "note")
classificationDF$classification = factor(classificationDF$classification, levels = c("Pathogenic", "Likely pathogenic", "VOUS", "Not classified", "Likely benign", "Benign", "Artifact")) 


#"leftGene", "leftSplicing", "leftAA", "rightGene", "rightSplicing", "rightAA",
annoted <- left_join(predSV, classificationDF, by = c("leftChr", "leftPos", "leftStrand", "rightChr", "rightPos", "rightStrand", "SVtype")) %>% 
  select(sample, sumReads, indelSize, leftChr, leftPos, leftStrand, NumofLeftSoftClippedReads, 
         rightChr, rightPos, rightStrand, NumofRightSoftClippedReads, SVtype, 
         leftGene, leftSplicing, leftAA, rightGene, rightSplicing, rightAA, classification, popAF, note, consensusSequences, everything()) %>% 
  replace_na(list(classification = "Not classified")) %>% 
  mutate(consensusSequences = sub("^","https://genome.ucsc.edu/cgi-bin/hgBlat?type=BLAT%27s+guess&userSeq=", consensusSequences)) %>% 
  arrange(classification, desc(sumReads))

openxlsx::write.xlsx(annoted, file = args[3])