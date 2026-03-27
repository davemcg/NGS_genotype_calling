#use with the new version of the variant scramble db.
args <- commandArgs(trailingOnly=TRUE)
library(tidyverse)
library(readxl)
#args <- c("Z:/exome/BlueprintGenetics/scramble_anno/all.exome.scramble.txt", "Z:/exome/BlueprintGenetics/scramble_anno/all.exome.scramble.v1.xlsx", "Z:/resources/SCRAMBLEvariantClassification.GRCh38.xlsx", "Z:/exome/BlueprintGenetics/scramble_anno/db.test.xlsx")
scramble_file <- args[1] 
output_xlsx_file <- args[2]
db_file <- args[3]
updated_db_file <- args[4]

scramble <- read_tsv(scramble_file, col_names = TRUE, na = c("NA", "", "None", "NONE", "."), col_types = cols(.default = col_character())) %>% 
  type_convert() %>% 
  filter(!is.na(Insertion)) %>% 
  separate(Insertion, c("chr", "pos"), sep = ":", remove = FALSE, convert = TRUE) %>% 
  mutate(pos = round(pos, -2)) %>% 
  unite("temp_Insertion", chr, pos, sep = ":", remove = TRUE) %>% 
  unite("temp_variantID", temp_Insertion, MEI_Family, Insertion_Direction, sep = "-", remove = FALSE ) %>% 
  select(-temp_Insertion)

scramble_family_count <- scramble %>%
  separate(sample, c("familyID", "individualID"), sep = "_", remove = FALSE) %>% 
  select(familyID) %>% 
  distinct() %>% 
  nrow()

scramble_count <- scramble %>% 
  separate(sample, c("familyID", "individualID"), sep = "_", remove = FALSE) %>% 
  select(temp_variantID, familyID) %>% 
  distinct() %>% 
  group_by(temp_variantID) %>% 
  summarise(CohortFreq = n()/scramble_family_count, AC = n(), AN = scramble_family_count) %>% 
  unite("NaltP/NtotalP", AC, AN, sep = "/", remove = TRUE) %>% 
  ungroup()
#74 families as shown by select(familyID) %>% distinct()

db_readme <- read_xlsx(db_file, sheet = "readme", na = c("NA", "", "None", "NONE", "."))
db <- read_xlsx(db_file, sheet = "Variant", na = c("NA", "", "None", "NONE", ".")) %>% 
  type_convert() %>%
  separate(Insertion, c("chr", "pos"), sep = ":", remove = FALSE, convert = TRUE) %>% 
  mutate(pos = round(pos, -2)) %>% 
  unite("temp_Insertion", chr, pos, sep = ":", remove = TRUE) %>% 
  unite("temp_variantID", temp_Insertion, MEI_Family, Insertion_Direction, sep = "-", remove = FALSE ) %>% 
  select(-temp_Insertion)
  
 # unite('variant', Insertion:Insertion_Direction, sep='-', remove = FALSE) 

db$autoClassification = factor(db$autoClassification, levels = c("Pathogenic", "Likely pathogenic", "VOUS", "VUS", "Not classified", "Likely benign", "Benign", "Artifact")) 
db$manualClassification = factor(db$manualClassification, levels = c("Pathogenic", "Likely pathogenic", "VOUS", "VUS", "Not classified", "Likely benign", "Benign", "Artifact")) 

scramble_analysis <- left_join(scramble, scramble_count, by = "temp_variantID") %>% 
  mutate(CohortFreq = pmax(cohortAF, CohortFreq, na.rm = TRUE)) %>% 
  mutate(autoClassification = case_when(CohortFreq > 0.2 ~ "Benign",
                                        CohortFreq > 0.06 ~ "Likely benign",
                                        TRUE ~ "Not classified")) %>% 
  select(-cohortAF)

db_update <- rbind(db, scramble_analysis) %>% distinct(temp_variantID, .keep_all = TRUE)

openxlsx::write.xlsx(scramble_analysis, file = output_xlsx_file)

openxlsx::write.xlsx(list("Variant" = db_update, "readme" = db_readme), file = updated_db_file, firstRow = TRUE)





  
  





