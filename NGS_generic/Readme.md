# Creating GVCF file(s) for each NGS sample

Based on Snakemake

This is a template for processing *generic* NGS files you have. As there are several types of input you may have (cram/bam/fastq) and an *infinite* number of ways they can be named and organized the Snakefile will have to be tweaked for your usage. 

The output of this Snakemake pipeline are some QC stats (fastqc, GATK metrics) and a GVCF file for each sample. 

The current Snakefile is wrapped by [Snakemake.wrapper.sh](https://github.com/davemcg/NGS_genotype_calling/blob/master/NGS_generic/Snakemake.wrapper.sh), which will run the Snakefile pipeline on [biowulf2](hpc.nih.gov).

Each file will be processed by chromosome (chr1 through X,Y) individually. chrMT and the contigs are merged together and processed as one. This makes this pipeline *fairly* performant. A exome can be processed in hours and a WGS in less than 48 hours. 

A major speed increase could be realized by subdividing the chromosome in (smaller pieces)[https://gatkforums.broadinstitute.org/gatk/discussion/10215/intervals-and-interval-lists] - the Broad breaks WGS into something like 500 pieces. But this would dramatically increase the complexity of the Snakemake pipeline and with [biowulf2](hpc.nih.gov) having limits, may be a bit stressful on the scheduler. 

It would be reasonable to split the autosomal chromosomes into, say, 2-4 pieces each. Something to do in the future. 



