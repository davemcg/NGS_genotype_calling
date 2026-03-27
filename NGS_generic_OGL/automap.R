args <- commandArgs(trailingOnly=TRUE)

# args <- c("Z:/exome/blueprint/AutoMap/G4V9_1_BP93806/G4V9_1_BP93806.HomRegions.tsv",
#           "Z:/exome/blueprint/AutoMap/test.annotSV.tsv", "Z:/resources/OGLpanelGeneDxORcandidate.xlsx", "Z:/exome/blueprint/AutoMap/test.annotSV.output.tsv")

automap_file <- args[1]
annotsv_file <- args[2]
geneCategory_file <- args[3]
annotated_file <- args[4]

library(tidyverse)
library(readxl)

automap <- read_tsv(automap_file, col_names = TRUE, na = c("NA", "", "None", "."), col_types = cols(.default = col_character())) %>%
  type_convert() %>% unite("key", `#Chr`:End, sep = "-", remove = FALSE) %>% mutate(key = sub("^chr", "", key))

annotsv <- read_tsv(annotsv_file, col_names = TRUE, na = c("NA", "", "None", "."), col_types = cols(.default = col_character())) %>%
  type_convert() %>% filter(Annotation_mode == "full") %>%
  mutate(SV_start = SV_start - 1) %>% 
  unite("key", SV_chrom:SV_end, sep = "-", remove = FALSE) %>% 
  select(key, Gene_name, Gene_count, OMIM_ID, OMIM_morbid, OMIM_morbid_candidate, DDD_HI_percent)

eyeGeneList <- read_xlsx(geneCategory_file, sheet = "analysis", na = c("NA", "", "None", ".")) %>% select(gene, panel_class) %>% pull(gene)

selectEyeGene <- function(x){
  x = as.character(x)
  geneNames <- as.list(strsplit(x, ";"))[[1]] 
  if (length(geneNames) > 0) {
    eyeGene <- purrr::keep(geneNames, geneNames %in% eyeGeneList) 
    if (length(eyeGene) == 0) {
      return(NA)
    } else if (length(eyeGene) == 1) { return(eyeGene) }
    else { return(paste(eyeGene, collapse = ";")) }
  } else {
    return(NA)
  }
}

annotsv$eyeGene <- sapply(1:nrow(annotsv), function(x) {selectEyeGene(annotsv[x, "Gene_name"])})

annotatedOutput <- left_join(automap,annotsv, by = "key") %>% select(`#Chr`:Percentage_homozygosity, eyeGene, Gene_name, Gene_count,OMIM_ID, OMIM_morbid, OMIM_morbid_candidate, DDD_HI_percent)

write_tsv(annotatedOutput, file.path('.', annotated_file))