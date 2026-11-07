#!/bin/bash
#SBATCH --gres=lscratch:100
#SBATCH --cpus-per-task=12
#SBATCH --mem=24g
#SBATCH --partition=norm,quick
#SBATCH --time=4:0:0

#seems not using lscratch, time = 1-1.5h
set -euo pipefail

###############################################################################
# alignment_to_fastq_by_rg.sh
#
# Convert BAM or CRAM to per-read-group FASTQ.gz files and create a mapping CSV.
#
# Usage:
#
#   BAM:
#     ./alignment_to_fastq_by_rg.sh \
#         input.bam \
#         bam2fastq \
#         mapping.csv
#
#   CRAM:
#     ./alignment_to_fastq_by_rg.sh \
#         input.cram \
#         bam2fastq \
#         mapping.csv \
#         old_reference.fa
#
# Output FASTQs:
#   <RG_ID>_R1_001.fastq.gz
#   <RG_ID>_R2_001.fastq.gz
#
# Mapping CSV columns:
#   sample,fastq,read_group
#
# The read_group field contains the complete original @RG line with
# literal tabs converted to \t, ready for:
#
#   bwa-mem2 mem -R "$read_group"
###############################################################################

input=$1
fastq_dir=$2
mapping_csv=$3
old_ref=${4:-}
module load GATK/4.7.0.0 samtools/1.23
tmpdir="/lscratch/${SLURM_JOB_ID:-$$}"

mkdir -p "$fastq_dir"

###############################################################################
# Check input type
###############################################################################

case "$input" in

    *.bam)
        gatk_ref_args=()
        samtools_ref_args=()
        ;;

    *.cram)

        pg_line=$(samtools view -H $input | grep '^@PG')
        case "$pg_line" in
          *"NCBI/GRCh38Decoy/genome.fa"*)
            old_ref="/data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa"
            ;;
          *"GATK"*"Homo_sapiens_assembly38.fasta"*)
            old_ref="/data/OGL/resources/genomes/GATK_resource_bundle/Homo_sapiens_assembly38.fasta"
            ;;
          *)
            echo "$old_ref"
            ;;
        esac
        gatk_ref_args=(-R "$old_ref")
        samtools_ref_args=(-T "$old_ref")
        ;;

    *)
        echo "ERROR: input must end in .bam or .cram" >&2
        exit 1
        ;;

esac


###############################################################################
# BAM/CRAM -> FASTQ.gz, one pair per read group
###############################################################################

echo "Converting:"
echo "  $input"
echo

gatk \
    --java-options "-Xmx20g -Djava.io.tmpdir=$tmpdir" \
    SamToFastq \
    -I "$input" \
    "${gatk_ref_args[@]}" \
    --RG_TAG ID \
    --OUTPUT_PER_RG true \
    --COMPRESS_OUTPUTS_PER_RG true \
    --OUTPUT_DIR "$fastq_dir" \
    --TMP_DIR "$tmpdir"

###############################################################################
# Rename:
#
#   RG_1.fastq.gz -> RG_R1_001.fastq.gz
#   RG_2.fastq.gz -> RG_R2_001.fastq.gz
###############################################################################

echo "Renaming FASTQ files..."

for fq in "$fastq_dir"/*_1.fastq.gz; do
    [[ -e "$fq" ]] || continue

    mv -- \
        "$fq" \
        "${fq%_1.fastq.gz}_R1_001.fastq.gz"
done

for fq in "$fastq_dir"/*_2.fastq.gz; do
    [[ -e "$fq" ]] || continue

    mv -- \
        "$fq" \
        "${fq%_2.fastq.gz}_R2_001.fastq.gz"
done

###############################################################################
# Extract original header
###############################################################################

header="$tmpdir/header.sam"

samtools view \
    "${samtools_ref_args[@]}" \
    -H \
    "$input" \
    > "$header"

###############################################################################
# Create mapping CSV
###############################################################################

#echo 'sample,fastq,read_group' > "$mapping_csv"

awk -F'\t' '$1=="@RG"' "$header" |
while IFS= read -r rg_line; do

    ###########################################################################
    # RG ID
    ###########################################################################

    rgid=$(
        printf '%s\n' "$rg_line" |
        awk -F'\t' '
        {
            for (i=1; i<=NF; i++) {
                if ($i ~ /^ID:/) {
                    sub(/^ID:/, "", $i)
                    print $i
                    exit
                }
            }
        }'
    )

    ###########################################################################
    # Sample name
    ###########################################################################

    sample=$(
        printf '%s\n' "$rg_line" |
        awk -F'\t' '
        {
            for (i=1; i<=NF; i++) {
                if ($i ~ /^SM:/) {
                    sub(/^SM:/, "", $i)
                    print $i
                    exit
                }
            }
        }'
    )

    if [[ -z "$rgid" ]]; then
        echo "ERROR: @RG entry does not contain ID:" >&2
        printf '%s\n' "$rg_line" >&2
        exit 1
    fi

    if [[ -z "$sample" ]]; then
        echo "ERROR: @RG $rgid does not contain SM" >&2
        exit 1
    fi

    ###########################################################################
    # Convert literal TABs in @RG line to \t for bwa-mem2 -R
    ###########################################################################

    bwa_rg=${rg_line//$'\t'/\\\\t}

    r1="${rgid}_R1_001.fastq.gz"
    r2="${rgid}_R2_001.fastq.gz"

    for fq in "$r1" "$r2"; do

        if [[ ! -f "$fastq_dir/$fq" ]]; then
            echo "WARNING: FASTQ not found: $fastq_dir/$fq" >&2
            continue
        fi

        printf '%s,%s,%s\n' \
            "$sample" \
            "$fq" \
            "$bwa_rg" \
            >> "$mapping_csv"

    done

done

###############################################################################
# Summary
###############################################################################

echo
echo "FASTQ files:"
find "$fastq_dir" \
    -maxdepth 1 \
    -type f \
    -name '*.fastq.gz' \
    -printf '%f\n' |
sort

echo
echo "Mapping CSV:"
cat "$mapping_csv"

echo
echo "Done."