
args <- commandArgs(trailingOnly=TRUE)
#When testing, comment out line above and use the line below.
 # args <- c("Z:/NextSeqAnalysis/221118/scramble_anno/G05256.forR.txt",
 #           "Y:/resources/SCRAMBLEvariantClassification.xlsx", "Z:/OGL_NGS/variant_prioritization/data/OGLv1_panel_DxORcandidate.tsv",
 #           "W:/ddl_nisc_custom_capture/042020/scramble/107910Y.annoted.xlsx")

#R_version: 'R/3.6.3'

library(tidyverse)
library(readxl)

annovar_file <- args[1]
scrambledb_file <- args[2]
geneCategory_file <- args[3]
HGMDtranscript_file <- args[4]
sampleName <- args[5]
output_file <- args[6]
output_xlsx_file <- args[7]

annovar <- read_tsv(annovar_file, col_names = TRUE, na = c("NA", "", "None", "."), col_types = cols(.default = col_character())) %>% 
  mutate(sample = sampleName) %>% 
  mutate(Intronic = gsub("0>-", "", Intronic)) %>% 
  unite("variant", Insertion:Insertion_Direction, sep = "-", remove = FALSE) %>% 
  separate(Insertion, c("chr", "pos"), sep = ":", remove = FALSE, convert = TRUE) %>% 
  mutate(pos = round(pos, -2)) %>% 
  unite("temp_Insertion", chr, pos, sep = ":", remove = TRUE) %>% 
  unite("temp_variantID", temp_Insertion, MEI_Family, Insertion_Direction, sep = "-", remove = FALSE ) %>% 
  select(-temp_Insertion)
 
# Do not use type_convert() because it converts chr:pos.

classificationDF <- read_xlsx(scrambledb_file, sheet = "Variant", na = c("NA", "", "None", ".")) %>% 
  select("temp_variantID", "autoClassification", "manualClassification", "note", "CohortFreq", "NaltP/NtotalP") %>% 
  distinct(temp_variantID, .keep_all = TRUE)

#panelGene <- read_tsv(args[3], col_names = TRUE, col_types = cols(.default = col_character())) %>% select(gene, panel_class)
panelGene <- read_xlsx(geneCategory_file, sheet = "analysis", na = c("NA", "", "None", ".")) %>% select(gene, panel_class)

HGMD <- read_tsv(HGMDtranscript_file, col_names = TRUE, na = c("NA", "", "None", "."), col_types = cols(.default = col_character())) 
hgmdNM <- dplyr::pull(HGMD, name)

mapHGMD <- function(x){
  x = as.character(x)
  annovarAA <- as.list(strsplit(x, ",|;"))[[1]] 
  if (length(annovarAA) > 1) {
    hgmdVar <- purrr::keep(annovarAA, str_extract_all(annovarAA, "NM_[[:digit:]]+") %in% hgmdNM) 
    if (length(hgmdVar) == 0) {
      return(x)
    } else if (length(hgmdVar) == 1) { return(hgmdVar) }
    else { return(paste(hgmdVar, collapse = ",")) }
  } else {
    return(x)
  }
}

annovar$AA <- sapply(1:nrow(annovar), function(x) {mapHGMD(annovar[x, "AA"])})
annovar$Intronic <- sapply(1:nrow(annovar), function(x) {mapHGMD(annovar[x, "Intronic"])})

#annovar <- annovar %>% unite("refgenewithver", Intronic, AA, sep = ",", remove = TRUE, na.rm = TRUE)

annoted <- left_join(annovar, classificationDF, by = c("temp_variantID")) %>% 
  replace_na(list(autoClassification = "Not classified", manualClassification = "Not classified")) 

annoted$autoClassification = factor(annoted$autoClassification, levels = c("Pathogenic", "Likely pathogenic", "VOUS", "VUS", "Not classified", "Likely benign", "Benign", "Artifact")) 
annoted$manualClassification = factor(annoted$manualClassification, levels = c("Pathogenic", "Likely pathogenic", "VOUS", "VUS", "Not classified", "Likely benign", "Benign", "Artifact")) 

chromosome_name <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM",
         "1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT")

annoted1 <- left_join(annoted, panelGene, by = c("Gene" = "gene")) %>%
  replace_na(list(panel_class = "Other")) %>% 
  mutate(eyeGene = case_when(panel_class == "Dx" ~ 9,
                             panel_class == "Candidate-High" ~ 8,
                             panel_class == "Candidate-refutedDxGene" ~ 7,
                             panel_class == "Candidate" ~ 6,
                             panel_class == "Dx-modifier-rare" ~ 3,
                             panel_class == "Dx-modifier-common" ~ 2,
                             TRUE ~ 0)) %>% 
  mutate(Clipped_Reads_In_Cluster = as.integer(Clipped_Reads_In_Cluster)) %>% 
  select(variant, Insertion, temp_variantID, MEI_Family, Insertion_Direction, Clipped_Reads_In_Cluster, Alignment_Score, 
         Alignment_Percent_Length, Alignment_Percent_Identity, Clipped_Sequence, Clipped_Side, Start_In_MEI, Stop_In_MEI, 
         polyA_Position, polyA_Seq, polyA_SupportingReads, TSD, TSD_length, panel_class, eyeGene, Func_refGene, Gene, Intronic, AA, autoClassification, manualClassification, sample, "note", "CohortFreq", "NaltP/NtotalP") %>% 
  separate(Insertion,
           c('temp_chromosome','temp_position'),
           sep = ':', remove = FALSE, convert = TRUE) %>% 
  filter(temp_chromosome %in% chromosome_name) %>% 
  arrange(desc(eyeGene), manualClassification, desc(Clipped_Reads_In_Cluster)) %>% 
  select(-temp_chromosome, -temp_position) 


#removed filter 6/30/2022, Func_refGene %in% c("splicing", "exonic", "UTR5", "UTR3", "upstream")
# [1] "downstream"          "intronic"            "exonic"              "intergenic"          "ncRNA_exonic"       
#[6] "ncRNA_intronic"      "ncRNA_splicing"      "splicing"            "upstream"            "UTR3"               
#[11] "upstream;downstream" "UTR5" 

if (dim(annoted1)[1] == 0) {
  print("###no scramble mei candidate after filtering with scramble db###")
  annoted1 <- annoted1 %>% add_row(note = "no scramble mei candidate after filtering with scramble db")
} 

write_tsv(annoted1, file.path('.',  output_file), na=".")

annoted2 <- filter(annoted1, !autoClassification %in% c("Likely benign", "Benign"),  !manualClassification %in% c("Likely benign", "Benign", "Artifact"), !Func_refGene %in% c("intergenic"))
if (dim(annoted2)[1] == 0) {
  print("###no scramble mei candidate after filtering with scramble db###")
  annoted2 <- annoted2 %>% add_row(note = "no scramble mei candidate after filtering with scramble db")
} 

openxlsx::write.xlsx(annoted2, file = output_xlsx_file)