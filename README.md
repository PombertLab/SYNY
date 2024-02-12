## <b>Synopsis</b>

The SYNY pipeline investigates gene colinearity (synteny) between genomes by reconstructing clusters from conserved pairs of protein-coding genes identified from [DIAMOND](https://github.com/bbuchfink/diamond) homology searches.

## <b>Table of contents</b>
* [Introduction](#Introduction)
* [Requirements](#Requirements)
  * [Downloading SYNY from GitHub](#downloading-SYNY-from-github)
  * [Installing dependencies](#installing-dependencies)
* [Using SYNY](#Using-SYNY)
  * [Command line options](#Command-line-options)
  * [A step by step example](#A-step-by-step-example)
* [References](#References)

## <b>Introduction</b>
#### <b>What is synteny?</b>
Synteny is the measure of similarity in gene organization between two species. Through time, genomes can be reorganized by large-scale events such as recombination, but also smaller events like -- but not exhaustively -- duplication, translocation, deletion, or inversion.

#### <b>Why use synteny?</b>
Synteny inferences can be used to:
- Improve genome annotation (by helping to identify mispredicted and/or missing genes)
- Perform evolutionary distance analyses (<i>e.g.</i> species sharing similar/identical genome reorganization events are likely more closely related than species that do not share these events).

## <b>Requirements</b>
- [DIAMOND](https://github.com/bbuchfink/diamond)
- [Perl5](https://www.perl.org/)
- [PerlIO::gzip](https://metacpan.org/pod/PerlIO::gzip)

##### <b>Optional</b>
- [Circos](https://circos.ca/)

#### <b>Downloading SYNY from GitHub</b>
```Bash
git clone https://github.com/PombertLab/SYNY.git
cd SYNY
export PATH=$PATH:$(pwd)
```

#### <b>Installing dependencies</b>
##### To install PerlIO::gzip:
```Bash
sudo cpan
  cpan[1]> install PerlIO::gzip
```

##### To install DIAMOND:
```Bash
version=v2.1.9     ## Replace with desired DIAMOND version
DIR=/opt/diamond   ## Replace with desired installation directory
mkdir -p $DIR

curl \
  -L https://github.com/bbuchfink/diamond/releases/download/$version/diamond-linux64.tar.gz \
  -o $DIR/diamond-linux64.tar.gz

tar -zxvf diamond-linux64.tar.gz --directory $DIR
rm diamond-linux64.tar.gz
export PATH=$PATH:$DIR
```

## <b>Using SYNY</b>
### Command line options
The SYNY pipeline can be run with [run_syny.pl](https://github.com/PombertLab/SYNY/blob/main/run_syny.pl), a master script that:
1. Extracts protein sequences from provided annotation files (currently supported: .gbf, .gbff, .prot, .embl).
2. Performs round-robin [DIAMOND](https://github.com/bbuchfink/diamond) BLASTP homology searches.
3. Generates gene pairs and reconstruct gene clusters.
4. Identifies sample-wide conserved genes, as well as missing or unique proteins for each species.
5. Generates configuration files/templates for plotting with [Circos](https://circos.ca/).

SYNY can be run from the master script as follows:<br>
```Bash
run_SYNY.pl \
  -a *.gbff \
  -g 5 \
  -o SYNY
```
Options for run_SYNY.pl are:
```
OPTIONS (MAIN):
-a (--annot)    Annotation files (Supported files: gbff, gff, embl)
-e (--evalue)   BLAST evalue cutoff [Default = 1e-10]
-g (--gaps)     Allowable number of gaps between pairs [Default = 0]
-o (--outdir)   Output directory [Default = SYNY]
-p (--prot)     Protein files # Generated automatically if GenBank files

OPTIONS (PLOTS): ##### Requires Circos - http://circos.ca/ #####
-r (--ref)      Genome to use as reference (defaults to first one
                alphabetically if none provided)
-u (--unit)     Size unit (Kb or Mb) [Default: Mb]
-c (--circos)   Generate Circos plots; currently buggy
                # works if run independently on configuration files
                # generated, e.g.: circos --conf concatenated.conf
-custom         Use custom color palette (20 colors) from Lemieux et al.:
                # https://pubmed.ncbi.nlm.nih.gov/31492891/
```
The output directory will be structured as follows: 
```Bash
drwxr-xr-x. 2 julian julian 4.0K May 10 18:18 ANNOTATIONS
drwxr-xr-x. 2 julian julian 4.0K May 10 18:18 CIRCOS
drwxr-xr-x. 2 julian julian 4.0K May 10 18:18 CONSERVED
drwxr-xr-x. 3 julian julian 4.0K May 10 18:18 DIAMOND
drwxr-xr-x. 2 julian julian 4.0K May 10 18:18 GENOME
drwxr-xr-x. 3 julian julian 4.0K May 10 19:29 LISTS
drwxr-xr-x. 2 julian julian 4.0K May 10 19:29 PROT_SEQ
drwxr-xr-x. 3 julian julian 4.0K May 10 18:18 SHARED
drwxr-xr-x. 3 julian julian 4.0K May 10 18:18 SYNTENY
```

The contents of the subdirectories are:
- ANNOTATIONS:
	- Tab-delimited lists of locus_tags and their products
- CIRCOS:
	- Configuration files for Circos plots
- CONSERVED:
	- Lists of proteins that are conserved within the samples (.conserved)
	- Summaries of all homology results and conservation within the samples (.conserved_summary)
	- Lists of unique proteins and their locations (.unique)
- DIAMOND:
	- DB
		- BLASTP databases
	- Round-robin BLASTP results (.diamond.6)
- GENOME:
	- FASTA files from the investigated genomes
- LISTS:
	- Lists of protein coding genes with location details (.list)
- PROT_SEQ:
	- Protein sequences for each species (.faa)
- SHARED:
	- Lists of all proteins and their top homologs (if any) in other species (.shared.tsv)
	- Lists of proteins that are unique to each species (.uniques.tsv)
- SYNTENY:
	- Cluster summary (clusters_summary.tsv)
	- Tab-delimited cluster summary table (clusters_summary_table.tsv)
	- Subdirectory per specified gap allowance (gap_#)
		- CLUSTERS:
			- Round-robin reconstructed syntetic clusters for each species
		- PAIRS:
			- Round-robin identified gene pairs for each species

### <b>A step by step example</b>

Below is a quick example describing how to compare two genomes from <i>Cryptococcus neoformans</i> var. <i>neoformans</i> JEC21 and <i>Cryptococcus gattii</i> WM276 using annotation data available in public databases.

##### Downloading annotation data from GenBank (NCBI):
```bash
DATA=~/DATA      ## Replace by desired annotation data directory
mkdir -p $DATA

# Cryptococcus neoformans JEC21
curl \
  -L https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/091/045/GCF_000091045.1_ASM9104v1/GCF_000091045.1_ASM9104v1_genomic.gbff.gz \
  -o $DATA/JEC21.gbff.gz

# Cryptococcus gattii WM276
curl \
  -L https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/185/945/GCF_000185945.1_ASM18594v1/GCF_000185945.1_ASM18594v1_genomic.gbff.gz \
  -o $DATA/WM276.gbff.gz
```

##### Running SYNY:
```Bash
SYNY=~/SYNY      ## Replace by desired SYNY output directory

run_syny.pl \
  -a $DATA/*.gbff.gz \
  -g 0 1 5 \
  -e 1e-10 \
  -r JEC21 \
  -o $SYNY
```

##### Plotting comparisons with Circos:
```Bash
CIRCOS=~/CIRCOS  ## Replace by desired Circos output directory

circos \
  -conf $SYNY/CIRCOS/concatenated/concatenated.conf \
  -outputdir $CIRCOS \
  -outputfile WM276_vs_JEC21.png
```

##### Example of clusters identified with SYNY
```Bash
head -n 29 $SYNY/SYNTENY/clusters_summary.tsv

##### JEC21_vs_WM276; Gap = 0 #####
  Total number of protein clusters:	5,760
  Longest:	64
  Shortest:	2
  Average cluster size:	7
  Median cluster size:	5
  N50:	10
  N75:	6
  N90:	4

##### JEC21_vs_WM276; Gap = 1 #####
  Total number of protein clusters:	5,921
  Longest:	204
  Shortest:	2
  Average cluster size:	26
  Median cluster size:	15
  N50:	53
  N75:	26
  N90:	13

##### JEC21_vs_WM276; Gap = 5 #####
  Total number of protein clusters:	5,956
  Longest:	434
  Shortest:	2
  Average cluster size:	112
  Median cluster size:	58
  N50:	237
  N75:	163
  N90:	66
```

```Bash
head -n 22 $SYNY/SYNTENY/gap_0/CLUSTERS/JEC21_vs_WM276.clusters

### Cluster 0001; CNA00020 to CNA00070; CGB_B0020C to CGB_B0070W; size = 6 ###
CNA00020	-	CGB_B0020C	-
CNA00030	+	CGB_B0030W	+
CNA00040	-	CGB_B0040C	-
CNA00050	+	CGB_B0050W	+
CNA00060	-	CGB_B0060C	-
CNA00070	+	CGB_B0070W	+
### Cluster 0002; CNA00080 to CNA00150; CGB_B0080W to CGB_B0150C; size = 8 ###
CNA00080	+	CGB_B0080W	+
CNA00090	+	CGB_B0090W	+
CNA00100	-	CGB_B0100C	-
CNA00110	+	CGB_B0110W	+
CNA00120	-	CGB_B0120C	-
CNA00130	+	CGB_B0130W	+
CNA00140	-	CGB_B0140C	-
CNA00150	-	CGB_B0150C	-
### Cluster 0003; CNA00180 to CNA00220; CGB_B0180C to CGB_B0220W; size = 5 ###
CNA00180	-	CGB_B0180C	-
CNA00190	-	CGB_B0190C	-
CNA00200	+	CGB_B0200W	+
CNA00210	-	CGB_B0210C	-
CNA00220	+	CGB_B0220W	+
```


##### Example of an image generated with Circos and SYNY (using defaults):
<p align="left">
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/WM276_vs_JEC21.png">
</p>


## <b>References</b>
[Sensitive protein alignments at tree-of-life scale using DIAMOND](https://www.nature.com/articles/s41592-021-01101-x). <b>Buchfink B, Reuter K, Drost HG.</b> Nature Methods 18, 366–368 (2021). doi:10.1038/s41592-021-01101-x

[Circos: an information aesthetic for comparative genomics](https://pubmed.ncbi.nlm.nih.gov/19541911/). <b >Krzywinski M, Schein J, Birol I, Connors J, Gascoyne R, Horsman D, Jones SJ, Marra MA.</b> Genome Res. 2009 Sep;19(9):1639-45. doi:10.1101/gr.092759.109