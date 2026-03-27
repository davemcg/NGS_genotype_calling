args <- commandArgs(trailingOnly=TRUE)

library(tidyverse)
library(readxl)

metadata_file <- args[1]
metadata_tsv_file <- args[2]

#prepare to write for multiple samples of the same target
metadata <- read_xlsx(metadata_file, sheet = "Primer", na = c("NA", "", "None", "NONE", ".")) %>% 
  select(Sample,Target,FwdPrimer,RevPrimer)

write_tsv(metadata, file.path('.',  metadata_tsv_file), col_names = FALSE, na=".")
