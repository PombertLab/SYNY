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

## Encephalitozoon intestinalis ATCC 50506 telomere-to-telomere (T2T) genome
outfile=${DATA}/Ei50506.gbff.gz
printf "\nDownloading Encephalitozoon intestinalis ATCC 50506 as ${outfile}\n\n"
curl \
  -L ${BASEURL}/GCA/024/399/295/GCA_024399295.1_ASM2439929v1/GCA_024399295.1_ASM2439929v1_genomic.gbff.gz \
  -o $outfile

## Encephalitozoon hellem ATCC 50604 telomere-to-telomere (T2T) genome
outfile=${DATA}/Eh50604.gbff.gz
printf "\nDownloading Encephalitozoon hellem ATCC 50604 as ${outfile}\n\n"
curl \
  -L ${BASEURL}/GCA/024/399/255/GCA_024399255.1_ASM2439925v1/GCA_024399255.1_ASM2439925v1_genomic.gbff.gz \
  -o $outfile

## Encephalitozoon cuniculi ATCC 50602 telomere-to-telomere (T2T) genome
outfile=${DATA}/Ec50602.gbff.gz
printf "\nDownloading Encephalitozoon cuniculi ATCC 50602 as ${outfile}\n\n"
curl \
  -L ${BASEURL}/GCA/027/571/585/GCA_027571585.1_ASM2757158v1/GCA_027571585.1_ASM2757158v1_genomic.gbff.gz \
  -o $outfile

## Ordospora colligata OC4
outfile=${DATA}/OcOC4.gbff.gz
printf "\nDownloading Ordospora colligata OC4 as ${outfile}\n\n"
curl \
  -L ${BASEURL}/GCF/000/803/265/GCF_000803265.1_ASM80326v1/GCF_000803265.1_ASM80326v1_genomic.gbff.gz \
  -o $outfile

## Ordospora pajunii FI-F-10
outfile=${DATA}/OpFIF10.gbff.gz
printf "\nDownloading Ordospora pajunii FI-F-10 as ${outfile}\n\n"
curl \
  -L ${BASEURL}/GCF/021/821/965/GCF_021821965.1_FI-F-10_v._1/GCF_021821965.1_FI-F-10_v._1_genomic.gbff.gz \
  -o $outfile
