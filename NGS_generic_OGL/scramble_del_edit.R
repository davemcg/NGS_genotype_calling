args <- commandArgs(trailingOnly=TRUE)

library(tidyverse)
library(readxl)

scramble_file <- args[1]
scrambledb_file <- args[2]
output_file <- args[3]

scramble <- read_tsv(scramble_file, col_names = TRUE, na = c("NA", "", "None", "."), col_types = cols(.default = col_character())) %>%
  type_convert() %>% 
  mutate(temp_SV_start = round(SV_start, -3), temp_SV_end = round(SV_end, -3)) %>%
  unite("variant", SV_chrom, temp_SV_start, temp_SV_end, SV_type, sep = "-", remove = FALSE) %>% 
  select(-temp_SV_start, -temp_SV_end)

scrambbleDB <- read_tsv(scrambledb_file, col_names = TRUE, na = c("NA", "", "None", "."), col_types = cols(.default = col_character())) %>% 
  type_convert()

annoted <- left_join(scramble, scrambbleDB, by = c("variant")) %>% 
  select(-variant)

write_tsv(annoted, file.path('.',  output_file), na=".")
