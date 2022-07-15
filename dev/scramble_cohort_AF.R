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
  unite('variant', Insertion:Insertion_Direction, sep='-', remove = FALSE) 

scramble_family_count <- scramble %>%
  separate(sample, c("familyID", "individualID"), sep = "_", remove = FALSE) %>% 
  select(familyID) %>% 
  distinct() %>% 
  nrow()

scramble_count <- scramble %>% 
  separate(sample, c("familyID", "individualID"), sep = "_", remove = FALSE) %>% 
  select(variant, familyID) %>% 
  distinct() %>% 
  group_by(variant) %>% 
  summarise(newCohortAF = n()/scramble_family_count )
#74 families as shown by select(familyID) %>% distinct()

db_readme <- read_xlsx(db_file, sheet = "readme", na = c("NA", "", "None", "NONE", "."))
db <- read_xlsx(db_file, sheet = "Variant", na = c("NA", "", "None", "NONE", ".")) %>% 
  type_convert() #%>% 
 # unite('variant', Insertion:Insertion_Direction, sep='-', remove = FALSE) 

db$autoClassification = factor(db$autoClassification, levels = c("Pathogenic", "Likely pathogenic", "VOUS", "Not classified", "Likely benign", "Benign", "Artifact")) 
db$manualClassification = factor(db$manualClassification, levels = c("Pathogenic", "Likely pathogenic", "VOUS", "Not classified", "Likely benign", "Benign", "Artifact")) 

scramble_analysis <- left_join(scramble, scramble_count, by = "variant") %>% 
  mutate(cohortAF = pmax(cohortAF, newCohortAF, na.rm = TRUE)) %>% 
  mutate(autoClassification = case_when(cohortAF > 0.2 ~ "Benign",
                                        cohortAF > 0.06 ~ "Likely benign",
                                        TRUE ~ "Not classified")) %>% 
  select(-newCohortAF)

db_cohortAF <- left_join(db, scramble_count, by = "variant") %>% 
  mutate(cohortAF = pmax(cohortAF, newCohortAF, na.rm = TRUE)) %>% 
  mutate(autoClassification = case_when(cohortAF > 0.2 ~ "Benign",
                                    cohortAF > 0.06 ~ "Likely benign",
                                    TRUE ~ "Not classified")) %>% 
  select(-newCohortAF)

manual_classified <-  filter(db_cohortAF, !is.na(manualClassification))
auto_classified <-  filter(db_cohortAF, is.na(manualClassification))

db_scramble <- bind_rows(auto_classified, scramble_analysis) %>% 
  group_by(variant) %>% 
  slice(which.max(cohortAF))

db_update <- bind_rows(manual_classified, db_scramble) %>% 
  distinct(variant, .keep_all = TRUE)

db_for_annotation <- db_update %>% 
  select("Insertion", "MEI_Family", "Insertion_Direction", "autoClassification", "manualClassification", "cohortAF", "note") 

scramble_partial <- scramble_analysis %>% select(variant, Insertion:AA, eyeGene, sample)
scramble_out <- left_join(scramble_partial, db_for_annotation, by = c("Insertion", "MEI_Family", "Insertion_Direction")) %>% 
  select(variant, Insertion, MEI_Family, Insertion_Direction, Clipped_Reads_In_Cluster, Alignment_Score, 
         Alignment_Percent_Length, Alignment_Percent_Identity, Clipped_Sequence, Clipped_Side, Start_In_MEI, Stop_In_MEI, 
         polyA_Position, polyA_Seq, polyA_SupportingReads, TSD, TSD_length, panel_class, eyeGene, Func_refGene, Gene, Intronic, AA, autoClassification, manualClassification, cohortAF, sample, note)

openxlsx::write.xlsx(scramble_out, file = output_xlsx_file)

openxlsx::write.xlsx(list("variant" = db_update, "readme" = db_readme), file = updated_db_file, firstRow = TRUE)





  
  





