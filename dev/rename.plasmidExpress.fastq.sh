#rename file 001_RGR015-002_reads.fastq.gz to RGR015-002_reads.fastq.gz
cd fastq
ls > ../rename.sh
cd ..
cut -d "_" -f 2- rename.sh | paste -d " " rename.sh - | sed 's/^/mv /' > fastq/rename.sh
cd fastq
bash rename.sh
rm rename.sh ../rename.sh
