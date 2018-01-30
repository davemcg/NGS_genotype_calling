# Creating GVCF file(s) for each NGS sample

Based on Snakemake

This is a template for processing *generic* NGS files you have. As there are several types of input you may have (cram/bam/fastq) and an *infinite* number of ways they can be named and organized the Snakefile will have to be tweaked for your usage. 

The output of this Snakemake pipeline are some QC stats (fastqc, GATK metrics) and a GVCF file for each sample. 

The current Snakefile is wrapped by [Snakemake.wrapper.sh](https://github.com/davemcg/NGS_genotype_calling/blob/master/NGS_generic/Snakemake.wrapper.sh), which will run the Snakefile pipeline on [biowulf2](hpc.nih.gov).

Each file will be processed by chromosome (chr1 through X,Y) individually. chrMT and the contigs are merged together and processed as one. This makes this pipeline *fairly* performant. A exome can be processed in hours and a WGS in less than 48 hours. 

A major speed increase could be realized by subdividing the chromosome in [smaller pieces](https://gatkforums.broadinstitute.org/gatk/discussion/10215/intervals-and-interval-lists) - the Broad breaks WGS into something like 500 pieces. But this would dramatically increase the complexity of the Snakemake pipeline and may be a bit stressful on the scheduler. 

It would be reasonable to split the autosomal chromosomes into, say, 2-4 pieces each. Something to do in the future. 

# Stuff to tweak for actual use
SAMPLES are currently identified by Snakemake's glob_wildcards function which scans the directory that the [Snakemake.wrapper.sh](https://github.com/davemcg/NGS_genotype_calling/blob/master/NGS_generic/Snakemake.wrapper.sh) is running in for files with the wildcard constraints given in the `wildcard_constraints:` rule in the [Snakefile](https://github.com/davemcg/NGS_genotype_calling/blob/master/NGS_generic/Snakefile). 

Depending on the naming scheme of your input files, you may have to give the Snakefile the actual file names in the Snakefile. Or feed it a file identified in the [config.yaml](https://github.com/davemcg/NGS_genotype_calling/blob/master/NGS_generic/config.yaml). 

It also takes fastq.gz as input. If you have a bam (very common), then you can look at `rule split_original_bam_by_rg:` and `rule align:` in the [EGA_EGAD00001002656_NGS_reanalyze Snakefile](https://github.com/davemcg/EGA_EGAD00001002656_NGS_reanalyze/blob/master/Snakefile) for inspiration on how to modify this Snakefile. 


