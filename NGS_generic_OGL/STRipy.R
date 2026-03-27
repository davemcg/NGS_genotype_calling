args <- commandArgs(trailingOnly=TRUE)
# args <- c("Z:/validation/genome-STR/STRipy/all.STRipy.tsv", 
#           "Z:/validation/genome-STR/STRipy/all.STRipy.xlsx")
tsv_file <- args[1]
excel_output_file <- args[2]

library(tidyverse)
library(readxl)

df <- read_tsv(tsv_file, col_names = TRUE, na = c("NA", "", ".", "None"), col_types = cols(.default = col_character())) %>%
  mutate(Allele1_CI= sub("^", " ", Allele1_CI), 
         Allele2_CI= sub("^", " ", Allele2_CI)) %>%
  type_convert() 
openxlsx::write.xlsx(list("STRipy" = df), file = excel_output_file, firstRow = TRUE, firstCol = TRUE)

