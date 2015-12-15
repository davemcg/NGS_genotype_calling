#!/bin/bash

#Prepends a single line onto a file

[ $# -eq 0 ] && { echo -e "\nPrepends a header to a file given.\
	 					   \nUsage: prepend.sh header_text file_to_prepend.\
						   \nWill overwrite file.\n"; exit 1; }

header=$1
file=$2

echo $header | cat - $file > ~/tmp/out && mv ~/tmp/out $file
