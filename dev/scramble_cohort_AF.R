#use with earlier version of the variant scramble db.
args <- commandArgs(trailingOnly=TRUE)
library(tidyverse)
library(readxl)
#args <- c("Z:/exome/BlueprintGenetics/scramble_anno/all.exome.scramble.txt", "Z:/exome/BlueprintGenetics/scramble_anno/all.exome.scramble.v2.xlsx", "Z:/resources/SCRAMBLEvariantClassification.GRCh38.xlsx", "Z:/resources/scramble/scrabmle.mei.db.xlsx")
scramble_file <- args[1] 
output_xlsx_file <- args[2]
db_file <- args[3]
updated_db_file <- args[4]

scramble <- read_tsv(scramble_file, col_names = TRUE, na = c("NA", "", "None", "NONE", "."), col_types = cols(.default = col_character())) %>% 
  type_convert() %>% 
  filter(!is.na(Insertion)) %>%
  unite('variant', Insertion:Insertion_Direction, sep='-', remove = FALSE) %>% 
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
  unite('variant', Insertion:Insertion_Direction, sep='-', remove = FALSE) %>% 
  separate(Insertion, c("chr", "pos"), sep = ":", remove = FALSE, convert = TRUE) %>% 
  mutate(pos = round(pos, -2)) %>% 
  unite("temp_Insertion", chr, pos, sep = ":", remove = TRUE) %>% 
  unite("temp_variantID", temp_Insertion, MEI_Family, Insertion_Direction, sep = "-", remove = FALSE ) %>% 
  select(-temp_Insertion)  

db$autoClassification = factor(db$autoClassification, levels = c("Pathogenic", "Likely pathogenic", "VOUS", "VUS", "Not classified", "Likely benign", "Benign", "Artifact")) 
db$manualClassification = factor(db$manualClassification, levels = c("Pathogenic", "Likely pathogenic", "VOUS", "VUS", "Not classified", "Likely benign", "Benign", "Artifact")) 

scramble_analysis <- left_join(scramble, scramble_count, by = "temp_variantID") %>% 
  mutate(CohortFreq = pmax(cohortAF, CohortFreq, na.rm = TRUE)) %>% 
  mutate(autoClassification = case_when(CohortFreq > 0.2 ~ "Benign",
                                        CohortFreq > 0.06 ~ "Likely benign",
                                        TRUE ~ "Not classified")) %>% 
  select(-cohortAF)

db_cohortAF <- left_join(db, scramble_count, by = "temp_variantID") %>% 
  mutate(CohortFreq = pmax(cohortAF, CohortFreq, na.rm = TRUE)) %>% 
  mutate(autoClassification = case_when(CohortFreq > 0.2 ~ "Benign",
                                        CohortFreq > 0.06 ~ "Likely benign",
                                        TRUE ~ "Not classified")) %>% 
  select(-cohortAF)

manual_classified <-  filter(db_cohortAF, !is.na(manualClassification))
auto_classified <-  filter(db_cohortAF, is.na(manualClassification))

db_scramble <- bind_rows(auto_classified, scramble_analysis) %>% 
  group_by(temp_variantID) %>% 
  slice(which.max(CohortFreq))

db_update <- bind_rows(manual_classified, db_scramble) %>% 
  distinct(temp_variantID, .keep_all = TRUE)

db_for_annotation <- db_update %>% 
  select(temp_variantID, "autoClassification", "manualClassification", "CohortFreq", "note") 

scramble_partial <- scramble_analysis %>% select(variant, temp_variantID, Insertion:AA, eyeGene, sample)
scramble_out <- left_join(scramble_partial, db_for_annotation, by = c("temp_variantID")) %>% 
  select(variant, Insertion, MEI_Family, Insertion_Direction, Clipped_Reads_In_Cluster, Alignment_Score, 
         Alignment_Percent_Length, Alignment_Percent_Identity, Clipped_Sequence, Clipped_Side, Start_In_MEI, Stop_In_MEI, 
         polyA_Position, polyA_Seq, polyA_SupportingReads, TSD, TSD_length, panel_class, eyeGene, Func_refGene, Gene, Intronic, AA, autoClassification, manualClassification, CohortFreq, sample, note)

openxlsx::write.xlsx(scramble_out, file = output_xlsx_file)

db_update_final <- db_update %>% select(-temp_variantID, -classification) %>% select(variant, everything())
openxlsx::write.xlsx(list("Variant" = db_update_final, "readme" = db_readme), file = updated_db_file, firstRow = TRUE)





  
  





