#use this if all of the samples have one underscore in filename.
for fastq1 in fastq/*.gz; do
	filename=$(basename $fastq1)
	header=$(zcat $fastq1 | head -1)
	id=$(echo $header | cut -d: -f 3,4 | sed 's/\:/\./g')
	sm=$(echo $filename | cut -d_ -f 1-2 | sed 's/\-/\_/g')
	echo "$sm,$filename,@RG\\\tID:$id"_"$sm\\\tSM:$sm\\\tLB:$lib"_"$sm\\\tPL:ILLUMINA" >> metadata_file.csv
	done
	
#for fastq1 in LPicturatus/*.gz; do filename=$(basename $fastq1); header=$(zcat $fastq1 | head -1); id=$(echo $header | cut -d: -f 3,4 | sed 's/\:/\./g'); sm="LPicturatus"; lib="nisc"; echo "$sm,$filename,@RG\\\tID:$id"_"$sm\\\tSM:$sm\\\tLB:$lib"_"$sm\\\tPL:ILLUMINA" >> metadata_file.csv; done