args <- commandArgs(trailingOnly=TRUE)
# args <- c("Z:/resources/bcmlocus.xlsx", "D1888-01", "Z:/projects/bcm_long_reads/annotation/D1888-01.avinput.hg38_multianno.txt", "Z:/projects/bcm_long_reads/annotation/D1888-01.bcm.test.tsv",
#           "Z:/resources/OGLpanelGeneDxORcandidate.xlsx", "rearranged.tsv", "filtered.tsv", "108976P", "filtered.xlsx", "0.5", "W:/ddl_nisc_custom_capture/042020/CoNVaDING/CNV_hiSens/108976P.b37.aligned.only.best.score.shortlist.txt")

tsv_file <- args[1]
excel_output_file <- args[2]

library(tidyverse)
library(readxl)

bcm <- read_tsv(tsv_file, col_names = TRUE, na = c("NA", "", ".", "None"), col_types = cols(.default = col_character())) %>%
  type_convert()

plof_path <- filter(bcm, grepl("splicing", Func.refGeneWithVer) 
                        | grepl("^frameshift|stop|start", ExonicFunc.refGeneWithVer) 
                        | grepl("path", ACMG_Class, ignore.case = TRUE))

openxlsx::write.xlsx(list("pLoF" = plof_path, "bcm" = bcm), file = excel_output_file, firstRow = TRUE, firstCol = TRUE)

