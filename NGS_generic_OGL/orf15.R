args <- commandArgs(trailingOnly=TRUE)
# args <- c("Z:/resources/bcmlocus.xlsx", "D1888-01", "Z:/projects/bcm_long_reads/annotation/D1888-01.avinput.hg38_multianno.txt", "Z:/projects/bcm_long_reads/annotation/D1888-01.bcm.test.tsv",
#           "Z:/resources/OGLpanelGeneDxORcandidate.xlsx", "rearranged.tsv", "filtered.tsv", "108976P", "filtered.xlsx", "0.5", "W:/ddl_nisc_custom_capture/042020/CoNVaDING/CNV_hiSens/108976P.b37.aligned.only.best.score.shortlist.txt")

tsv_file <- args[1]
excel_output_file <- args[2]

library(tidyverse)
library(readxl)

annovar <- read_tsv(tsv_file, col_names = TRUE, na = c("NA", "", ".", "None"), col_types = cols(.default = col_character())) %>%
  mutate(INFO=recode(INFO, "P" = "pileup", "F" = "full-alignment")) %>%
  separate(GT_FIELDS, c("GT","GQ","DP","AD","AF","PEAD"), sep=":", remove = FALSE) %>% 
  unite("refgenewithver", GeneDetail.refGeneWithVer, AAChange.refGeneWithVer, sep = ",", remove = TRUE, na.rm = TRUE) %>% 
  type_convert() %>% 
  select(Sample, CHROM:GT_FIELDS, Note, AD, PEAD, `Gene.refGeneWithVer`, Func.refGeneWithVer, ExonicFunc.refGeneWithVer, refgenewithver)
plof <- filter(annovar, grepl("splicing", Func.refGeneWithVer) | grepl("^frameshift|stop|start", ExonicFunc.refGeneWithVer) )

openxlsx::write.xlsx(list("pLoF" = plof, "orf15" = annovar), file = excel_output_file, firstRow = TRUE, firstCol = TRUE)

