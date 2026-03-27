args <- commandArgs(trailingOnly=TRUE)

# args <- c("Z:/exome/blueprint/AutoMap/G4V9_1_BP93806/G4V9_1_BP93806.HomRegions.tsv",
#           "Z:/exome/blueprint/AutoMap/test.annotSV.tsv", "Z:/resources/OGLpanelGeneDxORcandidate.xlsx", "Z:/exome/blueprint/AutoMap/test.annotSV.output.tsv")

mutserv_file <- args[1]
mitomap_file <- args[2]
annotated_file <- args[3]

library(tidyverse)


mitomap <- read_tsv(mitomap_file, col_names = TRUE, na = c("NA", "", "None", "."), col_types = cols(.default = col_character())) %>%
  select(-CHROM) %>% 
  mutate(Pubmed = gsub(",", ";", Pubmed)) %>% 
  type_convert()

mutserv <- read_tsv(mutserv_file, col_names = TRUE, na = c("NA", "", "None", "."), col_types = cols(.default = col_character())) %>%
  type_convert() %>% 
  mutate(ID = sub(".markDup.bam", "", ID)) %>% 
  left_join(., mitomap, by = c("Pos" = "POS", "Ref" = "REF", "Variant" = "ALT") )

write_tsv(mutserv, file.path('.', annotated_file), na="")