args <- commandArgs(trailingOnly=TRUE)

# args <- c("Z:/exome/blueprint/AutoMap/G4V9_1_BP93806/G4V9_1_BP93806.HomRegions.tsv",
#           "Z:/exome/blueprint/AutoMap/test.annotSV.tsv", "Z:/resources/OGLpanelGeneDxORcandidate.xlsx", "Z:/exome/blueprint/AutoMap/test.annotSV.output.tsv")

jax_file <- args[1]
sample_name <- args[2]
edited_file <- args[3]

library(tidyverse)

jax <- read_tsv(jax_file, col_names = TRUE, na = c("NA", "", "None", "."), col_types = cols(.default = col_character())) %>%
  type_convert() %>% 
  rename(SVtype = `user#1`, CN = `user#2`) %>% 
  mutate(SV_length = SV_end - SV_start) %>% 
  mutate(Samples_ID = sample_name)

write_tsv(jax, file.path('.', edited_file), na="")
