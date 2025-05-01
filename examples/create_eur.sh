#!/bin/bash

# File: create_files.bash
# Date: 1 May 2025
# Author: T. Quinn Smith
# Purpose: Create VCF file from 1000 Genomes Project for an example.

# Print out commands while it executes.
set -eux

# The samples are in EUR_40.txt.
# 10 samples from IBS, GBR, FIN, and TSI in that order.

# ALL_chr6.vcf.gz contains all samples from the 1000 Genomes Project Chromosome 6.
# Get portion of p-arm and remove monomorphic sites.
tabix ALL_chr6.vcf.gz
bcftools view -r "chr6:1-40000000" -S EUR_40.txt ALL_chr6.vcf.gz | bcftools view -c 1 -c 1:nonmajor | gzip > EUR_40_chr6.vcf.gz

# Create target file.
wget https://raw.githubusercontent.com/gavinr/world-countries-centroids/refs/heads/master/dist/countries.csv

awk -F ',' '{
	R = 6371;
	PI = 3.1415926535;
	long = $1 * PI / 180;
	lati = $2 * PI / 180;
	x = R * cos(lati) * cos(long);
	y = R * cos(lati) * sin(long);
	if ($3 == "Spain" || $3 == "United Kingdom" || $3 == "Finland" || $3 == "Italy") {
		print $4"\t"x"\t"y;
	}
}' countries.csv > EUR_coords.tsv

sed "4q;d" EUR_coords.tsv | awk '{for (i=1;i<=10;i++) print $2"\t"$3}' > EUR_40_target.tsv
seq 1 10 | awk '{print "IBS"}' > EUR_40_labels.txt
sed "1q;d" EUR_coords.tsv | awk '{for (i=1;i<=10;i++) print $2"\t"$3}' >> EUR_40_target.tsv
seq 1 10 | awk '{print "GBR"}' >> EUR_40_labels.txt
sed "3q;d" EUR_coords.tsv | awk '{for (i=1;i<=10;i++) print $2"\t"$3}' >> EUR_40_target.tsv
seq 1 10 | awk '{print "FIN"}' >> EUR_40_labels.txt
sed "2q;d" EUR_coords.tsv | awk '{for (i=1;i<=10;i++) print $2"\t"$3}' >> EUR_40_target.tsv
seq 1 10 | awk '{print "TSI"}' >> EUR_40_labels.txt
