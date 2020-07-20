#the following files follow David's coverage R code. now downloaded to /data/OGL/resources/
#ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/GRCh37_mapping/gencode.v33lift37.metadata.HGNC.gz
#ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/GRCh37_mapping/gencode.v33lift37.basic.annotation.gtf.gz

library(data.table)
library(tidyverse)
library(readxl)

args <- commandArgs(trailingOnly=TRUE)
#When testing, comment out line above and use the line below.
#args <- c("Y:/NextSeqAnalysis/pipelineValid/139.thresholds.bed",
#         "Z:/OGL_NGS/variant_prioritization/data/OGLv1_panel_DxORcandidate.tsv", "test.xlsx")


readDepth <- fread(args[1], header=T) %>% mutate(targetSize = end - start) 
#panelGene <- read_tsv(args[2], col_names = TRUE, col_types = cols(.default = col_character())) %>% select(gene, panel_class)
panelGene <- read_xlsx(args[2], sheet = "analysis", na = c("NA", "", "None", ".")) %>% select(gene, panel_class)

RDtype <- left_join(readDepth, panelGene, by = c("gene")) %>% replace_na(list(panel_class="Other"))
RDtype$panel_class = factor(RDtype$panel_class, levels = c("Dx", "Candidate", "Other"))

totalRegion <- RDtype %>% group_by(panel_class) %>% summarise(totalRegionSum = sum(targetSize))
coverageTen <- RDtype %>% group_by(panel_class) %>% summarise(coverageTenSum = sum(coverageTen)) %>% select(coverageTenSum)
coverageTwenty <- RDtype %>% group_by(panel_class) %>% summarise(coverageTwentySum = sum(coverageTwenty)) %>% select(coverageTwentySum)
coverageThirty <- RDtype %>% group_by(panel_class) %>% summarise(coverageThirtySum = sum(coverageThirty)) %>% select(coverageThirtySum)

percentage <- cbind.data.frame(totalRegion, coverageTen, coverageTwenty, coverageThirty) %>% 
  mutate(percentTen = coverageTenSum/totalRegionSum, percentTwenty = coverageTwentySum/totalRegionSum, percentThirty = coverageThirtySum/totalRegionSum) 

lessTen <- RDtype %>% filter(coverageTen < targetSize) %>% arrange(panel_class)
lessTwenty <- RDtype %>% filter(coverageTwenty < targetSize) %>% arrange(panel_class)
lessThirty <- RDtype %>% filter(coverageThirty < targetSize) %>% arrange(panel_class)

openxlsx::write.xlsx(list("PercentCoverage" = percentage, "less10Xregion" = lessTen, "less20Xregion" = lessTwenty, "less30Xregion" = lessThirty), file = args[3])

