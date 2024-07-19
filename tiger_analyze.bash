#!/bin/bash
 
# This is the analysis script for Part 2 of the TIGeR Bioinformatician Interview Assignment
# Author : Viraj Muthye
# Date: 7/18/2024

# The automated workflow in Part 1 has identified one potential problematic VCF file 
# and has saved the names of files without errors in a file named "VCF_files_without_errors".
# This script loops through all the valif VCF files identified in Part 1 by the automated workflow

# usage : ./tiger_analyze.bash

# Prerequisites that need to be downloaded prior to running the workflow
# 	- vcftools (similar to Part 1)
#	- SnpEff & SnpSift
#	- Java

# Installation instructions for vcftools via GitHub
# git clone https://github.com/vcftools/vcftools.git
# cd vcftools/
# ./autogen.sh
# export PERL5LIB=/path/to/your/vcftools-directory/src/perl/
# ./configure prefix=$HOME
# make
# make install
# export PATH=/path/to/vcf/bin/folder/:$PATH *
# export PERL5LIB=vcftools/src/perl/ *
# These last two steps marked with an asterisk (*) are extremely important and the cause of several downstream issues if not done

# Installation for SnpEff & SnpSift
# Install this in the same directory you will run this script from **
# wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
# unzip snpEff_latest_core.zip
# This is all we need to do to install. We now need to download the human hg38 genome database.
# Download the Human Hg38 genome database : java -jar snpEff/snpEff.jar download GRCh38.86 

# Installation for Java
# Go to http://java.com and click on the Download button. https://www.java.com/en/download/linux_manual.jsp
# cd directory_path_name
# tar zxvf jre-8u73-linux-x64.tar.gz
# When the installation has completed, you will see the word Done
# More instructions here : https://www.java.com/en/download/help/linux_x64_install.html

# Output Files
# <FILENAME>.vcf.gz_Report.txt                         - The report file which includes the unique SNP calls and summary of multiple features
# <FILENAME>.vcf.gz_unique_snps.recode.vcf.gz          - VCF file with the unique SNPs
# <FILENAME>.vcf.gz_FS_histogram.png		- PNG file of a histogram of Phred-scaled p-value using Fisher's exact test
# <FILENAME>.vcf.gz_allele_frequency_histogram.png	- PNG file of a histogram of allele frequencies
# <FILENAME>.vcf.gz_unique_snps.recode.annotated.vcf_missense_variants.vcf.gz - VCF file of missense variants

while read p;do 
	
	echo "This is the report for the file : $p." > "$p"_Report.txt
	
	echo "These are the unique SNP calls identified in $p. Below this long list, I have added information on multiple features of this VCF file." >> "$p"_Report.txt
	echo "################################################" >> "$p"_Report.txt
	
	# sort the vcf file
	vcftools --vcf "$p"_FILTERED.recode.vcf --recode --recode-INFO-all --stdout | vcf-sort > "$p"_sorted_snps.vcf

	# get the unique positions of the SNPs 
	grep -v "^#" "$p"_sorted_snps.vcf | awk '{print $1"\t"$2}' | uniq -u > "$p"_unique_snps_positions.txt
	
	# filter the original file to keep the unique SNPs
	vcftools --vcf "$p"_sorted_snps.vcf --positions "$p"_unique_snps_positions.txt --recode --recode-INFO-all --out "$p"_unique_snps

	
	# add unique SNPs to report
	echo "Number of unique SNP calls in $p." >> "$p"_Report.txt
	vcftools --vcf "$p"_sorted_snps.vcf --positions "$p"_unique_snps_positions.txt --recode --recode-INFO-all --stdout | grep -v "#" | awk '{print $1"\t"$2"\t"$4"\t"$5}' | sort | uniq | wc -l >> "$p"_Report.txt
	echo " " >> "$p"_Report.txt

	vcftools --vcf "$p"_sorted_snps.vcf --positions "$p"_unique_snps_positions.txt --recode --recode-INFO-all --stdout | grep -v "#" | awk '{print $1"\t"$2"\t"$4"\t"$5}' >> "$p"_Report.txt
	echo "The VCF file with the uniqe SNP calls is $p_unique_snps.recode.vcf" >> "$p"_Report.txt 	
	echo " " >> "$p"_Report.txt
	              
	# Now, we move on to the features
	# Feauture 1
	echo "Feature 1: Chromosomal Distribution of SNPs" >> "$p"_Report.txt
	echo "################################################" >> "$p"_Report.txt
	echo "Distribution of SNPs on autosomes and sex chromosomes will impact our biological interpretation of the SNPs and also will also help us filter our data for downstream analyses. In many instances, we would want to remove SNPs on sex chromosomes, but several diseases are sex-linked and for such analyses we would want to include the SNPs on chromosomes X and Y. In fact, a vast proportion of GWAS studies exclude SNPs on the X chromosome, with some of the issues cited being "complications in genotype calling, imputation, selection of test statistics" (PMID: 24408308)." >> "$p"_Report.txt
	
	# extract the information from column 1 which contains the chrosome number, sort it, and add the counts per chromosome to the report
	awk '{print $1}' "$p"_unique_snps.recode.vcf | grep -v "#" | sort | uniq -c | sort -k1nr >> "$p"_Report.txt
	echo " " >> "$p"_Report.txt
	
	# Identifying variants that cause missense mutations that affect coding regions
	echo "Of particular interest are missense mutations. These are mutations that alter the genetic code in a way that it induces a change in the amino acid in that particular position. This is just one of the different ways we can filter variants by using the different effects in the ANN field of the annotations." >> "$p"_Report.txt
	
	java -jar snpEff/snpEff.jar GRCh38.86 "$p"_unique_snps.recode.vcf > "$p"_unique_snps.recode.annotated.vcf
	java -jar snpEff/SnpSift.jar filter "ANN[*].EFFECT has 'missense_variant'" "$p"_unique_snps.recode.annotated.vcf > "$p"_unique_snps.recode.annotated.vcf_missense_variants.vcf
	rm "$p"_unique_snps.recode.annotated.vcf
	echo "Number of SNPs that are mapped to missense variants:" >> "$p"_Report.txt
	grep -v "#" "$p"_unique_snps.recode.annotated.vcf_missense_variants.vcf | awk '{print $1"\t"$2"\t"$4"\t"$5}' | sort | uniq | wc -l  >> "$p"_Report.txt
	echo " " >> "$p"_Report.txt

	echo "Additionally, it would be informative to know where the SNPs are on each chromosome to know more about the biological relevance of these regions. Here, we analyzed the SNP density in a bin size of 5,000 bases on just one chromosome - chromosome X. We present the top 25 most dense bins in this report." >> "$p"_Report.txt
	# extract the data from the VCF file for the X chromosome
	vcftools --vcf "$p"_unique_snps.recode.vcf --chr chrX --recode --recode-INFO-all --out "$p"_chrX_snps

	# calculate the SNP density in bins of 5,000 bases on the X chromosome and save it to a temporary file
	vcftools --vcf "$p"_chrX_snps.recode.vcf --SNPdensity 5000 --out "$p"_chrX_snps_snpd
	head -1 "$p"_chrX_snps_snpd.snpden >> "$p"_Report.txt

	# remove empty bins, sort by the size of the bin, and then export the data to the report
	cat "$p"_chrX_snps_snpd.snpden | awk '$3>0{print $1"\t"$2"\t"$3"\t"$4}' | sort -k3nr | grep -v "CHROM" | head -25 >> "$p"_Report.txt
	echo " " >> "$p"_Report.txt


	# Feature 2
	echo "Feature 2: Detecting strand bias of SNPs using Phred-scaled p-value using Fisher's exact test" >> "$p"_Report.txt
	echo "################################################" >> "$p"_Report.txt	
	echo "I chose this feature as an example for analyzing the quality of variant calls. Strand bias is a form of sequencing bias where one DNA strand is favored over the other, resulting in potential incorrect evaluation of the amount of evidence observed for one allele. One way to identify this in our data is through the FisherStrand annotation. This method uses Fisher's Exact Test to determine if there is strand bias between forward and reverse strands for the reference or alternate allele. The output of this method is a Phred-scaled p-value, which we can use to analyze the data. More bias is usually indicative of false positive calls." >> "$p"_Report.txt
	echo " " >> "$p"_Report.txt

	# extract the p-values from the vcf file 
	vcftools --vcf "$p"_unique_snps.recode.vcf --get-INFO FS --out tmpFS

	# plot the distribution of p-values and save it as a PNG image
	Rscript plot_FS.R
	
	# identify average p-value, maximim p=value, and the number of SNPs with a p-value more than 0
	average=$(awk '{print $5}' tmpFS.INFO | grep -v "FS" | awk '{total += $1; count++ } END {print total/count}')
	maximum=$(awk '{print $5}' tmpFS.INFO | grep -v "FS" | sort -k1nr | head -1)
	above=$(awk '{print $5}' tmpFS.INFO | grep -v "FS" | grep -v "0" | sort | uniq | wc -l)
	echo "The average p-value is $average, the maximum value is $maximum, and $above SNPs have a p-value more than 0." >> "$p"_Report.txt

	# delete the INFO file
	rm tmpFS.INFO

	# rename the image file
	mv allele_FS_histogram.png "$p"_FS_histogram.png
	echo " " >> "$p"_Report.txt    
	echo "Check out the distribution of the p-values in the file $p _FS_histogram.png" >> "$p"_Report.txt

	echo " " >> "$p"_Report.txt

	# Feauture 3 : Number of Biallelic and Multiallelic SNPs 
	echo "Feature 3: Number of biallelic and multiallelic SNPs" >> "$p"_Report.txt
	echo "################################################" >> "$p"_Report.txt  	  
	echo "This data contains both biallelic and multiallelic SNPs. A biallelic site is a specific locus in a genome that contains two observed alleles. A multiallelic site is a specific locus in a genome that contains three or more observed alleles. True instances of multiallelic sites are not very common unless you look at very large cohorts, so such instances are many times taken as a sign of a noisy region where artefacts could be likely. It is important to identify and we might want to filter these sites for downstream analysis. For instance, most GWAS studies will remove these multiallelic variants, however multiple multiallelic variants have been identified as disease relevant (PMID: 27535533 and PMID: 26670213)." >> "$p"_Report.txt
	echo " " >> "$p"_Report.txt

	# extract the allele frequencies of SNPs 
	vcftools --gzvcf "$p"_unique_snps.recode.vcf --freq2 --out "$p"_unique_snps.recode.vcf_freq
	echo "Number of bi-allelic SNPs" >> "$p"_Report.txt
	tail -n+2 "$p"_unique_snps.recode.vcf_freq.frq | awk '$3==2 {print $1"\t"$2}' | sort | uniq | wc -l >> "$p"_Report.txt
	echo "Number of multi-allelic SNPs" >> "$p"_Report.txt
	tail -n+2 "$p"_unique_snps.recode.vcf_freq.frq | awk '$3>2  {print $1"\t"$2}' | sort | uniq | wc -l >> "$p"_Report.txt
	echo " " >> "$p"_Report.txt   
	
	# Distribution of allele frequency of bilallelic SNPs
	echo "Allele frequency is an essential characteristic of SNP analysis as it tells us how prevalent a SNP is. SNP frequency is determined by the frequency of second most prevalent allele : the minor allele frequency. The frequency of an allele within a population and the risk for complex diseases associated with that allele are important characteristics to consider in downstream analyses." >> "$p"_Report.txt
	echo col1$'\t'col2$'\t'col3$'\t'col4$'\t'col5 > tmp.txt
	tail -n+2 "$p"_unique_snps.recode.vcf_freq.frq | awk '$3==2 {print $1"\t"$2"\t"$3"\t"$4"\t"$5}'>> tmp.txt	
	
	# this R script will plot the distribution of allele frequencies and save the output in a PNG file
	Rscript plot_AF.R

	# delete the INFO file
	rm tmp.txt

	# rename the PNG file
	mv allele_frequency_histogram.png "$p"_allele_frequency_histogram.png
	
	echo " " >> "$p"_Report.txt
	echo "The distrbution of p-values can be viewed in the file $p _allele_frequency_histogram.png" >> "$p"_Report.txt	

done < VCF_files_without_errors

for file in *.vcf;do gzip $file;done

# delete intermediate files and clean up the folder
rm *unique_snps_positions.txt
rm *.log
rm *chrX*
rm *sorted_snps.vcf.gz
rm *vcf_freq.frq
rm snpEff_genes.txt
rm snpEff_summary.html






