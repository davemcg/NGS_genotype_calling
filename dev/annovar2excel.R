args <- commandArgs(trailingOnly=TRUE)

# args <- c("Z:/exome/blueprint/AutoMap/G4V9_1_BP93806/G4V9_1_BP93806.HomRegions.tsv",
#           "Z:/exome/blueprint/AutoMap/test.annotSV.tsv", "Z:/resources/OGLpanelGeneDxORcandidate.xlsx", "Z:/exome/blueprint/AutoMap/test.annotSV.output.tsv")

annovar_file <- args[1]
output_file <- args[2]
#geneCategory_file <- args[3]
#annotated_file <- args[4]

library(tidyverse)
library(readxl)

annovar <- read_tsv(annovar_file, col_names = TRUE, na = c("NA", "", "None", "."), col_types = cols(.default = col_character())) %>%
  mutate(INFO = case_when( INFO == "P" ~ "pileup",
                           INFO == "F" ~ "full-alignment",
                           TRUE ~ INFO)) %>%
  type.convert() %>% 
  mutate(Note = ifelse(POS %in% c(38285414, 38298269, 38299739), "Homopolymer", "")) %>% 
  select(CHROM:GT_FIELDS, Note, `Func.refGeneWithVer`:`AAChange.refGeneWithVer`)

openxlsx::write.xlsx(list("orf15" = annovar), file = output_file, firstRow = TRUE, firstCol = FALSE)
