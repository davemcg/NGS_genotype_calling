library(ggplot2)
library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
#When testing, comment out line above and use the line below.

#args <- c("G02210.b37.aligned.only.normalized.coverage.txt", "file.pdf", "1.1", "0.8", "chr_gain_or_loss.txt")

CoNVaDING <- args[1] # normalized_coverage_from_CoNVaDING
chrRDfilename <- args[2]
highCutOff <- args[3]
lowCutOff <- args[4]
abnormalChr_filename <- args[5]


ChrMeanRD <- read.delim(CoNVaDING, sep = "\t", header = TRUE, colClasses = c("factor","integer","integer","factor","numeric","numeric","numeric","numeric",
                                  "numeric","numeric","numeric") ) %>%
  group_by(CHR) %>%
  summarize(avg_of_normalized_autosomal=mean(NORMALIZED_AUTOSOMAL))

ChrMeanRD$CHR <- factor(ChrMeanRD$CHR, levels=c(paste(1:22,sep=""),"X", "Y"))
ChrMeanRD = ChrMeanRD[order(ChrMeanRD$CHR),]

ChrMeanRD_plot <- ggplot(ChrMeanRD,aes(x=CHR,y=avg_of_normalized_autosomal))+ 
  #scale_color_discrete (breaks=c("1","2","3","4","5","6","7","8","9","10", "11", "12","13","14","15","16","17","19","20", "21","22","X"))+
  geom_hline(yintercept = 1, color = "red", linetype="dotted")+
  geom_point()+
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2), limits=c(0,2)) +
  theme_minimal() +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=10)) +
  labs(x = "Chromosome", y = "Chromosome Autosomal-Normalized Read Depths")

pdf(file = chrRDfilename,   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches

# Step 2: Create the plot with R code
plot(ChrMeanRD_plot)

# Step 3: Run dev.off() to create the file!
dev.off()


AbnormalChr <- ChrMeanRD %>% filter(CHR != "X") %>%
  filter(avg_of_normalized_autosomal > highCutOff | avg_of_normalized_autosomal < lowCutOff)
#if AbnormalChr has a row, then print, Can consider sex from ped file then print plot only when the condition is met.
if (dim(AbnormalChr)[1] != 0) {
  write_tsv(AbnormalChr, file.path('.', abnormalChr_filename))
}

