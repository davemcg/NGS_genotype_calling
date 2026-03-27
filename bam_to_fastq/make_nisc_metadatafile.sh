lib=$1
for fastq1 in fastq/*.fastq.gz; do
	filename=$(basename $fastq1)
	header=$(zcat $fastq1 | head -1 | sed 's/@//')
	id=$(echo $header | cut -d: -f 1,2,8 | sed 's/\:/\./g')
	sm=$(echo $filename | cut -d_ -f 1 | sed 's/\-/\_/g')
	echo "$sm,$filename,@RG\\\tID:$id"_"$sm\\\tSM:$sm\\\tLB:$lib"_"$sm\\\tPL:ILLUMINA" >> metadata_file.csv
	done
	
#C32GLACXX:8:2112:08299:15551 1:N:0:1
