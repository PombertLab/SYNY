#!/usr/bin/env bash

DATA=$1 ## Desired annotation data directory

## Readme
if [[ ( $@ == "--help") ||  $@ == "-h" || $# -le 0 ]]
then 
	printf "\nUsage: $0 output_directory/\n\n"
	exit 0
fi 

### Create dir
mkdir -p $DATA

##### Downloading data from NCBI ####
BASEURL=https://ftp.ncbi.nlm.nih.gov/genomes/all

### Arabidopsis thaliana TAIR10.1 - 119.1 Mbp
outfile=${DATA}/TAIR10.gbff.gz
printf "\nDownloading Arabidopsis thaliana TAIR10.1 (119.1 Mbp) as ${outfile}\n\n"
curl \
  -L ${BASEURL}/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.gbff.gz \
  -o $outfile

### Arabidopsis arenosa AARE701a - 149.7 Mbp
outfile=${DATA}/AARE701.gbff.gz
printf "\nDownloading Arabidopsis arenosa AARE701a (149.7 Mbp) as ${outfile}\n\n"
curl \
  -L ${BASEURL}/GCA/905/216/605/GCA_905216605.1_AARE701a/GCA_905216605.1_AARE701a_genomic.gbff.gz \
  -o $outfile
