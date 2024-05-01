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

### Arabidopsis arenosa AARE701a - 149.7 Mbp
outfile=${DATA}/AARE701.gbff.gz
printf "\nDownloading Arabidopsis arenosa AARE701a (149.7 Mbp) as ${outfile}\n\n"
curl \
  -L ${BASEURL}/GCA/905/216/605/GCA_905216605.1_AARE701a/GCA_905216605.1_AARE701a_genomic.gbff.gz \
  -o $outfile

### Arabidopsis suecica As9502 - 271.6 Mbp
outfile=${DATA}/As9502.gbff.gz
printf "\nDownloading Arabidopsis suecica As9502 (271.6 Mbp) as ${outfile}\n\n"
curl \
  -L ${BASEURL}/GCA/019/202/805/GCA_019202805.1_ASM1920280v1/GCA_019202805.1_ASM1920280v1_genomic.gbff.gz \
  -o $outfile

### Arabidopsis thaliana x Arabidopsis arenosa Allo738 - 268.6 Mbp
outfile=${DATA}/Allo738.gbff.gz
printf "\nDownloading Arabidopsis thaliana x Arabidopsis arenosa Allo738 (268.6 Mbp) as ${outfile}\n\n"
curl \
  -L ${BASEURL}/GCA/019/202/795/GCA_019202795.1_ASM1920279v1/GCA_019202795.1_ASM1920279v1_genomic.gbff.gz \
  -o $outfile

### Arabidopsis thaliana TAIR10.1 - 119.1 Mbp
outfile=${DATA}/TAIR10.gbff.gz
printf "\nDownloading Arabidopsis thaliana TAIR10.1 (119.1 Mbp) as ${outfile}\n\n"
curl \
  -L ${BASEURL}/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.gbff.gz \
  -o $outfile
