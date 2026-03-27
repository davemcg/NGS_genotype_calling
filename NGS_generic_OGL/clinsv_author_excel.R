args <- commandArgs(trailingOnly=TRUE)

# args <- c("Z:/resources/OGLpanelGeneDxORcandidate.xlsx",
#           "Z:/genome/RodYoung/clinSV/RY1.clinSV.RARE_PASS_GENE.annotated.tsv",
#           "Z:/development/genome/clinSV/D1596_1.RARE_PASS_GENE.eG.tsv",
#           "Z:/development/genome/clinSV/D1596_1.RARE_PASS_GENE.eG.filtered.xlsx")

geneCategory_file <- args[1]
heavy_file <- args[2]
edited_tsv_file <- args[3]
edited_xlsx_file <- args[4]

library(tidyverse)
library(readxl)

eyeGeneList <- read_xlsx(geneCategory_file, sheet = "analysis", na = c("NA", "", "None", ".")) %>% 
  select(gene) %>% distinct() %>% 
  pull(gene)

heavy <- read_xlsx(heavy_file, sheet = "Sheet1", 
                   col_names = c("Sample", "VariantID", "Filter", "isRare", "SU",
                                 "PAF_SU", "PE", "SR", "DRF", "DRA", "PAF_DRA", "PCSD",
                                 "GT", "MQBP", "isCNV", "IGV", "GOTO", "Genomic_Location",
                                 "SVtype", "SVlen", "Tool", "PopAF_MGRB", "Empty",
                                 "PopAF1k","GC", "CR", "MQ", "SEGD", "NumberOfGenes",	
                                 "Genes", "GeneFeature", "HPO", "PHEN"), na = c("NA", "", "None", ".")) 


selectEyeGene <- function(x){
  x = as.character(x)
  geneNames <- as.list(strsplit(x, ","))[[1]] 
  if (length(geneNames) > 0) {
    eyeGene <- purrr::keep(geneNames, geneNames %in% eyeGeneList) 
    if (length(eyeGene) == 0) {
      return(NA)
    } else if (length(eyeGene) == 1) { return(eyeGene) }
    else { return(paste(eyeGene, collapse = ",")) }
  } else {
    return(NA)
  }
}

heavy$eyeGene <- sapply(1:nrow(heavy), function(x) {selectEyeGene(heavy[x, "Genes"])})

heavy <- select(heavy, Sample:SEGD, eyeGene, everything()) %>% 
  arrange(eyeGene)
write_tsv(heavy, edited_tsv_file, na = "")

heavy_filtered <- heavy %>% 
  filter(grepl(",chr", Genomic_Location) | SVtype != "BND") 
openxlsx::write.xlsx(list("clinSV" = heavy_filtered), file = edited_xlsx_file, firstRow = TRUE, firstCol = TRUE)
