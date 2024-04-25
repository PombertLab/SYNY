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

### Cryptococcus neoformans strain JEC21
outfile=${DATA}/JEC21.gbff.gz
printf "\nDownloading Cryptococcus neoformans strain JEC21 as ${outfile}\n\n"
curl \
  -L ${BASEURL}/GCF/000/091/045/GCF_000091045.1_ASM9104v1/GCF_000091045.1_ASM9104v1_genomic.gbff.gz \
  -o $outfile

### Cryptococcus gattii strain WM276
outfile=${DATA}/WM276.gbff.gz
printf "\nDownloading Cryptococcus neoformans strain JEC21 as ${outfile}\n\n"
curl \
  -L ${BASEURL}/GCF/000/185/945/GCF_000185945.1_ASM18594v1/GCF_000185945.1_ASM18594v1_genomic.gbff.gz \
  -o $outfile

### Cryptococcus gattii VGV strain MF34
outfile=${DATA}/MF34.gbff.gz
printf "\nDownloading Cryptococcus gattii VGV strain MF34 as ${outfile}\n\n"
curl \
  -L ${BASEURL}/GCA/009/650/685/GCA_009650685.1_Cryp_gatt_MF34/GCA_009650685.1_Cryp_gatt_MF34_genomic.gbff.gz \
  -o $outfile

### Cryptococcus decagattii strain 7685027
outfile=${DATA}/D7685.gbff.gz
printf "\nDownloading Cryptococcus decagattii strain 7685027 as ${outfile}\n\n"
curl \
  -L ${BASEURL}/GCA/036/417/295/GCA_036417295.1_ASM3641729v1/GCA_036417295.1_ASM3641729v1_genomic.gbff.gz \
  -o $outfile

### Cryptococcus deuterogattii strain R265
outfile=${DATA}/R265.gbff.gz
printf "\nDownloading Cryptococcus deuterogattii strain R265 as ${outfile}\n\n"
curl \
  -L ${BASEURL}/GCF/002/954/075/GCF_002954075.1_C._deuterogattii_R265_chr/GCF_002954075.1_C._deuterogattii_R265_chr_genomic.gbff.gz \
  -o $outfile