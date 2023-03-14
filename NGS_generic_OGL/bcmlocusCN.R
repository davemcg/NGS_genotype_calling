args <- commandArgs(trailingOnly=TRUE)
#args <- c( "D1695_02", "Z:/genome/nisc23-1/coverage/mosdepth/D1695_02.md.mosdepth.summary.txt", "Z:/genome/nisc23-1/bcmlocus/mosdepth/D1695_02.md.regions.bed.gz", "Z:/genome/nisc23-1/bcmlocus/D1695_02.bcm.cn.tsv", "Z:/genome/nisc23-1/bcmlocus/D1695_02.bcm.cn.wide.tsv")

sampleName <- args[1]
chrCoverage_file <- args[2]
bcm_coverage_bed_file <- args[3]
bcmlocusCN_file <- args[4]
bcmlocusCN_wide_file <- args[5]
# CN calculation; rearrangedGemini_file <- args[3]

library(tidyverse)
library(readxl)
library(vroom)

chrCoverage <- read_tsv(chrCoverage_file, col_names = TRUE, na = c("NA", "", "None"), col_types = cols(.default = col_character())) %>% 
  slice(1:48) %>% 
  type_convert() %>% 
  filter(!grepl("_region", chrom)) %>% 
  rename(median_coverage = mean)
autosome <- chrCoverage %>% slice(1:22) %>% summarise(median_coverage = median(median_coverage)) %>% 
  mutate(target = "autosome", CN = 2) %>% 
  select(target, median_coverage, CN)
autosome_median_coverage <- autosome %>% pull(median_coverage)
sexChr <- filter(chrCoverage, chrom %in% c("chrX", "chrY") ) %>% 
  rename(target = chrom) %>% 
  select(target, median_coverage) %>% 
  mutate(CN = round(2*median_coverage/autosome_median_coverage, 1))

bcm_coverage_bed <- vroom(bcm_coverage_bed_file, col_names = FALSE) %>% 
  rename(chr = X1, start = X2, end = X3, target = X4, median_coverage = X5) %>% 
  separate(target, c("gene", "probeNo"), sep = "_", remove = FALSE) %>% 
  filter(!gene %in% c("TEX28", "OPN1MW2", "OPN1MW3")) #%>%
 # filter(!target %in% c("OPN1MW_2", "OPN1MW_3", "OPN1MW_4", "OPN1MW_5", "OPN1MW_6"))

bcm_locus <- filter(bcm_coverage_bed, gene %in% c("OPSIN-LCRup", "OPSIN-LCRcore", "OPSIN-LCRdn", "OPN1LW", "OPN1MW")) %>% 
  select(target, median_coverage) %>% 
  mutate(CN = round(2*median_coverage/autosome_median_coverage, 1))

x_non_bcm_locus <- filter(bcm_coverage_bed, !gene %in% c("OPSIN-LCRup", "OPSIN-LCRcore", "OPSIN-LCRdn", "OPN1LW", "OPN1MW")) %>% 
  group_by(gene) %>% 
  summarise(median_coverage = median(median_coverage)) %>% 
  rename(target = gene) %>% 
  mutate(CN = round(2*median_coverage/autosome_median_coverage, 1))

allTarget <- rbind(autosome, sexChr, x_non_bcm_locus, bcm_locus) %>% 
  mutate(target = factor(target, levels = c("autosome", "chrX", "chrY", "NAA10", "RENBP", "HCFC1", "TMEM187", "IRAK1", "MECP2", "OPSIN-LCRup_1", "OPSIN-LCRup_2", 
                                            "OPSIN-LCRup_3","OPSIN-LCRup_4","OPSIN-LCRcore", "OPSIN-LCRcore_1", "OPSIN-LCRdn_1", "OPSIN-LCRdn_2", "OPSIN-LCRdn_3", "OPSIN-LCRdn_4", "OPSIN-LCRdn_5",
                                            "OPN1LW_promoter","OPN1LW_1", "OPN1LW_In1","OPN1LW_2","OPN1LW_3","OPN1LW_4", "OPN1LW_5", "OPN1LW_6", "OPN1MW_promoter", "OPN1MW_1",
                                            "OPN1MW_2","OPN1MW_3","OPN1MW_4", "OPN1MW_5", "OPN1MW_6","TKTL1", "FLNA"))) %>% 
  arrange(target)


write_tsv(allTarget, file = bcmlocusCN_file)

allTargetWider <- allTarget %>% select(target, CN) %>% 
  pivot_wider(names_from = target, values_from = CN) %>% 
  mutate(sample = sampleName) %>% 
  select(sample, everything())

write_tsv(allTargetWider, file = bcmlocusCN_wide_file)



