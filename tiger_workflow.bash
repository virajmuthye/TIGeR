#!/bin/bash

# This is an automated workflow for Part 1 of the TIGeR Bioinformatician Interview Assignment
# Author : Viraj Muthye
# Date: 7/18/2024

# Prerequisites that need to be downloaded prior to running the workflow : vcftools
# git clone https://github.com/vcftools/vcftools.git
# cd vcftools/
# ./autogen.sh
# export PERL5LIB=/path/to/your/vcftools-directory/src/perl/
# ./configure prefix=$HOME
# make
# make install
# export PATH=/path/to/vcf/bin/folder/:$PATH

# Check if arguments have been provided, else print usage and quit
if [ $# -eq 0 ]; then
    echo "Usage: $0 arg1 arg2 ..."
    exit 1
fi

# If arguments have been provided, loop through all the arguments
for arg in "$@"
do

	# run vcf-validator and save the total number of errors in the variable num
	echo "File: $arg"
	num=$(vcf-validator -d -u $arg | grep "error" | awk '{print $1}')

	# filter the files based on the number of errors identified
	if [ $num -gt 0 ]; then

		# save the name of the vcf file with atleast one error in VCF_files_with_errors 
    		echo "$arg" >> VCF_files_with_errors 
	else

		# save the name of vcf file without any error in VCF_files_without_errors
		# then, remove all indels and keep just the SNPs, filter with a minimum Phred QUAL score of 20, filter with a minimum depth DP of 5 reads, and then output the SNPS to a new vcf file 
    		echo "$arg" >> VCF_files_without_errors
    		vcftools --gzvcf $arg --remove-indels --minQ 20 --minDP 10 --recode --recode-INFO-all --out "$arg"_FILTERED	
	fi

done















