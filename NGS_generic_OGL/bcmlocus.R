args <- commandArgs(trailingOnly=TRUE)
# args <- c("Z:/resources/bcmlocus.xlsx", "D1888-01", "Z:/projects/bcm_long_reads/annotation/D1888-01.avinput.hg38_multianno.txt", "Z:/projects/bcm_long_reads/annotation/D1888-01.bcm.test.tsv",
#           "Z:/resources/OGLpanelGeneDxORcandidate.xlsx", "rearranged.tsv", "filtered.tsv", "108976P", "filtered.xlsx", "0.5", "W:/ddl_nisc_custom_capture/042020/CoNVaDING/CNV_hiSens/108976P.b37.aligned.only.best.score.shortlist.txt")

bcmlocus_file <- args[1]
sampleName <- args[2]
annovar_file <- args[3]
annovar_edited_file <- args[4]
annovar_edited_excelFile <- args[5]
# CN calculation; rearrangedGemini_file <- args[3]

library(tidyverse)
library(readxl)

bcmlocus <- read_xlsx(bcmlocus_file, sheet = "bcmlocus", na = c("NA", "", "None", ".")) %>%
  select(ID,AAChange.refGeneWithVer,Annotation:ACMG_Class) %>% 
  separate(ID, c("caller", "variant"), sep = "_", remove = TRUE) %>% 
  separate(AAChange.refGeneWithVer, c("Gene", "Transcript", "Exon", "HGVSc", "HGVSp"), sep = ":", remove = TRUE, convert = TRUE) %>% 
  select(caller,Gene,HGVSc:ACMG_Class)

annovar <- read_tsv(annovar_file, col_names = TRUE, na = c("NA", "", "None", "."), col_types = cols(.default = col_character())) %>%
  type_convert() %>%
  separate_rows(Gene.refGeneWithVer, sep = ";") %>% 
  filter(!Gene.refGeneWithVer %in% c("OPN1MW2", "OPN1MW3")) %>% 
  separate_rows(AAChange.refGeneWithVer, sep = ",") %>% 
  filter(!grepl("OPN1MW2|OPN1MW3", AAChange.refGeneWithVer)) %>% 
  separate(AAChange.refGeneWithVer, c("Gene", "Transcript", "Exon", "HGVSc", "HGVSp"), sep = ":", remove = FALSE, convert = TRUE) %>% 
  mutate(Note = "") %>% 
  mutate(Sample = sampleName) %>% 
  filter(QUAL > 0) %>% 
  separate(ID, c("caller", "variant"), sep = "_", remove = FALSE)

annovar_edit <- left_join(annovar, bcmlocus, by=c("caller","Gene","HGVSc", "HGVSp")) %>%
  separate(GT_FIELDS, c("GT","DP","VAF","PS"), sep = ":", convert = TRUE) %>%
  separate(INFO, c("AO", "temp_DP"), sep = ";", remove = FALSE) %>% 
  mutate(AO = as.integer(sub("AO=", "", AO))) %>% 
  mutate(AO = ifelse(is.na(INFO), round(DP * VAF, 0), AO)) %>% 
  mutate(VAF=round(VAF, 2)) %>% 
  mutate(GT = ifelse(caller == "haplo", "", GT)) %>% 
  select(CHROM, POS, ID, REF, ALT, QUAL, Func.refGeneWithVer, Gene.refGeneWithVer, GeneDetail.refGeneWithVer, ExonicFunc.refGeneWithVer, AAChange.refGeneWithVer,
         HGVSp, Annotation, Function, ACMG_Class, Note, Sample,VAF,AO,DP,GT,PS) %>% 
  replace_na(list(Annotation = "", Function = "", ACMG_Class = ""))


# annovarE3 <- filter(annovar, Exon == "exon3")
# annovarE245 <- filter(annovar, Exon != "exon3")
# 
# annovarE245_edit <- left_join(annovarE245, bcmlocus, by=c("HGVSp")) %>%
#   select(CHROM, POS, ID, REF, ALT, QUAL, Func.refGeneWithVer, Gene.refGeneWithVer, GeneDetail.refGeneWithVer, ExonicFunc.refGeneWithVer, AAChange.refGeneWithVer,
#          HGVSp, Annotation, Function, ACMG_Class, Note, Sample, INFO, FORMAT, GT_FIELDS)
# 
# annovarE3_edit <- left_join(annovarE3, bcmlocus, by=c("HGVSc", "HGVSp")) %>% 
#   select(CHROM, POS, ID, REF, ALT, QUAL, Func.refGeneWithVer, Gene.refGeneWithVer, GeneDetail.refGeneWithVer, ExonicFunc.refGeneWithVer, AAChange.refGeneWithVer,
#          HGVSp, Annotation, Function, ACMG_Class, Note, Sample, INFO, FORMAT, GT_FIELDS)
# 
# annovar_edit <- rbind(annovarE245_edit, annovarE3_edit) %>%
#   arrange(POS) %>% 
#   replace_na(list(Annotation = "", Function = "", ACMG_Class = ""))
  
write_tsv(annovar_edit, file = annovar_edited_file, na=".")
openxlsx::write.xlsx(list("bcm" = annovar_edit), file = annovar_edited_excelFile, firstRow = TRUE, firstCol = TRUE)

