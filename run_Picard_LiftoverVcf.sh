#!/bin/bash

module load picard

input_vcf=$1
output_vcf_name=$2
java -Xmx8g -jar $PICARDJARPATH/picard.jar LiftoverVcf \
	I=$input_vcf \
	O=$output_vcf_name \
	CHAIN=
	
