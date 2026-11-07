#use this if all of the samples have one underscore in filename.
for fastq1 in fastq/*.gz; do
	filename=$(basename $fastq1)
	header=$(zcat $fastq1 | head -1)
	id=$(echo $header | cut -d: -f 3,4 | sed 's/\:/\./g')
	sm=$(echo $filename | cut -d_ -f 1_2 | sed 's/\-/\_/g' | sed 's/_/x/')
	nfields=$(awk -F: '{print NF}' <<< "$header")
	if [[ "$nfields" -lt 10 ]]; then
		pu="${id}_${sm}"
	elif [[ "$nfields" -ge 10 ]]; then
		pu=$(echo "$header" | cut -d: -f 3,4,10 | sed 's/:/./g')
	else
		echo "WARNING: unexpected FASTQ header format ($nfields fields): $header" >&2
		continue
	fi
	echo "$sm,$filename,@RG\\\tID:${id}_${sm}\\\tSM:$sm\\\tLB:${lib}_${sm}\\\tPU:$pu\\\tPL:ILLUMINA" \
		>> $metadata_file
	done

#for fastq1 in LPicturatus/*.gz; do filename=$(basename $fastq1); header=$(zcat $fastq1 | head -1); id=$(echo $header | cut -d: -f 3,4 | sed 's/\:/\./g'); sm="LPicturatus"; lib="nisc"; echo "$sm,$filename,@RG\\\tID:$id"_"$sm\\\tSM:$sm\\\tLB:$lib"_"$sm\\\tPL:ILLUMINA" >> metadata_file.csv; done
