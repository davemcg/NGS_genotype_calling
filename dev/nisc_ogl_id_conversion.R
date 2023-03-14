library(tidyverse)
library(readxl)
args <- commandArgs(trailingOnly=TRUE)

conversion_file <- args[1]
metadata_file <- args[2]
newMetadata_file <- args[3]

conversion <- read_tsv(conversion_file, col_names = TRUE, na = c("NA", "", "None", "NONE", ".", "FALSE", "False"), col_types = cols(.default = col_character())) %>%
  type_convert()

metadata <- read_csv(metadata_file, col_names = FALSE, na = c("NA", "", "None", "NONE", ".", "FALSE", "False"), col_types = cols(.default = col_character())) %>%
  type_convert() %>% 
  left_join(., conversion, by = c("X1" = "niscID")) 

write_csv(metadata, file = newMetadata_file, col_names = FALSE)

