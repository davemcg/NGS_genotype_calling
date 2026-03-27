#add each of three amplicon reads to the sample fastq files, so that PORECHOP will not exit with error of "no adaptor found". This appears to push PORECHOP to find the other real reads. These reads were removed from PORECHOP fastq output if output reads > 1.
#find the 18 kb RPGR amplicon in IGV
#find good LW and MW reads in run1 on IGV 
zgrep -n "b31e9787-cd66-4e79-a6a0-f63388dd6811" ../run1/fastq/D703x001.fastq.gz
zcat ../run1/fastq/D703x001.fastq.gz | sed -n '13073,13076p' | sed '1s/^@/@LW_EXAMPLE/' >> ~/git/NGS_genotype_calling/ont/bcmORF15.fastq
zgrep -n "64d58d32-2a65-4b77-911d-14b42cf7de45" ../run1/fastq/D703x001.fastq.gz
zcat ../run1/fastq/D703x001.fastq.gz | sed -n '16741,16744p' | sed '1s/^@/@MW_EXAMPLE/' >> ~/git/NGS_genotype_calling/ont/bcmORF15.fastq
zgrep -n "218ce401-751f-4b74-8df8-bc76d75e60d2" ../run1/fastq/D703x001.fastq.gz
zcat ../run1/fastq/D703x001.fastq.gz | sed -n '13885,13888p' | sed '1s/^@/@RPGR_EXAMPLE/' >> ~/git/NGS_genotype_calling/ont/bcmORF15.fastq
