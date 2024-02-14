## <b>Synopsis</b>

The SYNY pipeline investigates gene colinearity (synteny) between genomes by reconstructing clusters from conserved pairs of protein-coding genes identified from [DIAMOND](https://github.com/bbuchfink/diamond) homology searches.

## <b>Table of contents</b>
* [Introduction](#Introduction)
* [Requirements](#Requirements)
  * [Downloading SYNY from GitHub](#downloading-SYNY-from-github)
  * [Installing dependencies](#installing-dependencies)
* [Using SYNY](#Using-SYNY)
  * [Command line options](#Command-line-options)
  * [Step by step examples](#Step-by-step-examples)
    * [Example 1: Cryptococcus](#Example-1---Cryptococcus)
    * [Example 2: Encephalitozoon](#Example-2---Encephalitozoon)
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
	- Subdirectory per specified gap allowance (gap_#) containing:
		- CLUSTERS:
			- Round-robin reconstructed syntetic clusters for each species
		- PAIRS:
			- Round-robin identified gene pairs for each species

### <b>Step by step examples</b>

#### Example 1 - Cryptococcus
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
mkdir -p $CIRCOS

## Default orientation
circos \
  -conf $SYNY/CIRCOS/concatenated/concatenated.conf \
  -outputdir $CIRCOS \
  -outputfile WM276_vs_JEC21.png

## Inverted orientation
circos \
  -conf $SYNY/CIRCOS/concatenated/concatenated.inverted.conf \
  -outputdir $CIRCOS \
  -outputfile WM276_vs_JEC21.inverted.png
```

##### Example of clusters identified with SYNY
```Bash
head -n 32 $SYNY/SYNTENY/clusters_summary.tsv

##### JEC21_vs_WM276; Gap = 0 #####
  Total number of proteins in clusters:	5760
  # of clusters:	779
  Longest:	64
  Shortest:	2
  Average cluster size:	7
  Median cluster size:	5
  N50:	10
  N75:	6
  N90:	4

##### JEC21_vs_WM276; Gap = 1 #####
  Total number of proteins in clusters:	5921
  # of clusters:	228
  Longest:	204
  Shortest:	2
  Average cluster size:	26
  Median cluster size:	15
  N50:	53
  N75:	26
  N90:	13

##### JEC21_vs_WM276; Gap = 5 #####
  Total number of proteins in clusters:	5956
  # of clusters:	53
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

##### Example of an image generated with Circos and SYNY (using the inverted karyotype):
<p align="left">
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/WM276_vs_JEC21.inverted.png">
</p>


#### Example 2 - Encephalitozoon
Below is a quick example describing how to compare a total of three telomere-to-telomere (T2T) genomes from <i>Encephalitozoon</i> species <i>E. intestinalis</i> ATCC 50506, <i>E. hellem</i> ATCC 50451, and <i>E. cuniculi</i> ATCC 50602 using annotation data available in public databases.

##### Downloading annotation data from GenBank (NCBI):
```bash
DATA=~/ENCE
mkdir -p $DATA

## Encephalitozoon intestinalis ATCC 50506 telomere-to-telomere (T2T) genome
curl \
  -L https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/024/399/295/GCA_024399295.1_ASM2439929v1/GCA_024399295.1_ASM2439929v1_genomic.gbff.gz \
  -o $DATA/intestinalis_50506.gbff.gz

## Encephalitozoon hellem ATCC 50451 telomere-to-telomere (T2T) genome
curl \
  -L https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/029/215/505/GCA_029215505.1_ASM2921550v1/GCA_029215505.1_ASM2921550v1_genomic.gbff.gz \
  -o $DATA/hellem_50451.gbff.gz

## Encephalitozoon cuniculi ATCC 50602 telomere-to-telomere (T2T) genome
curl \
  -L https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/027/571/585/GCA_027571585.1_ASM2757158v1/GCA_027571585.1_ASM2757158v1_genomic.gbff.gz \
  -o $DATA/cuniculi_50602.gbff.gz

```

##### Running SYNY:
```Bash
SYNY=~/ENCE      ## Replace by desired SYNY output directory

run_syny.pl \
  -a $DATA/*.gbff.gz \
  -g 0 1 5 \
  -e 1e-10 \
  -r intestinalis_50506 \
  -o $SYNY
```

##### Plotting comparisons with Circos:
```Bash
CIRCOS=~/CIRCOS  ## Replace by desired Circos output directory
mkdir -p $CIRCOS

## Default orientation
circos \
  -conf $SYNY/CIRCOS/concatenated/concatenated.conf \
  -outputdir $CIRCOS \
  -outputfile encephalitozoon.png
```

##### Example of clusters identified with SYNY
```Bash
head -n 32 $SYNY/SYNTENY/clusters_summary.tsv

##### cuniculi_50602_vs_hellem_50451; Gap = 0 #####
  Total number of proteins in clusters: 1831
  # of clusters:        87
  Longest:      101
  Shortest:     2
  Average cluster size: 21
  Median cluster size:  14
  N50:  38
  N75:  22
  N90:  10

##### cuniculi_50602_vs_hellem_50451; Gap = 1 #####
  Total number of proteins in clusters: 1840
  # of clusters:        26
  Longest:      185
  Shortest:     2
  Average cluster size: 71
  Median cluster size:  49
  N50:  148
  N75:  86
  N90:  49

##### cuniculi_50602_vs_hellem_50451; Gap = 5 #####
  Total number of proteins in clusters: 1845
  # of clusters:        19
  Longest:      195
  Shortest:     2
  Average cluster size: 97
  Median cluster size:  92
  N50:  160
  N75:  147
  N90:  55
```

```Bash
head -n 22 $SYNY/SYNTENY/gap_0/CLUSTERS/intestinalis_50506_vs_cuniculi_50602.clusters

### Cluster 001; GPK93_01g00070 to GPK93_01g00480; J0A71_11g22950 to J0A71_11g23360; size = 41 ###
GPK93_01g00070  -       J0A71_11g22950  -
GPK93_01g00080  +       J0A71_11g22960  +
GPK93_01g00090  +       J0A71_11g22970  +
GPK93_01g00100  +       J0A71_11g22980  +
GPK93_01g00110  +       J0A71_11g22990  +
GPK93_01g00120  +       J0A71_11g23000  +
GPK93_01g00130  -       J0A71_11g23010  -
GPK93_01g00140  -       J0A71_11g23020  -
GPK93_01g00150  -       J0A71_11g23030  -
GPK93_01g00160  +       J0A71_11g23040  +
GPK93_01g00170  -       J0A71_11g23050  -
GPK93_01g00180  +       J0A71_11g23060  +
GPK93_01g00190  -       J0A71_11g23070  -
GPK93_01g00200  -       J0A71_11g23080  -
GPK93_01g00210  -       J0A71_11g23090  -
GPK93_01g00220  +       J0A71_11g23100  +
GPK93_01g00230  -       J0A71_11g23110  -
GPK93_01g00240  -       J0A71_11g23120  -
GPK93_01g00250  -       J0A71_11g23130  -
GPK93_01g00260  -       J0A71_11g23140  -
GPK93_01g00270  -       J0A71_11g23150  -
```


##### Example of an image generated with Circos and SYNY comparing a total of 3 genomes (using defaults):
<p align="left">
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/encephalitozoon.png">
</p>

## <b>References</b>
[Sensitive protein alignments at tree-of-life scale using DIAMOND](https://www.nature.com/articles/s41592-021-01101-x). <b>Buchfink B, Reuter K, Drost HG.</b> Nature Methods 18, 366–368 (2021). doi:10.1038/s41592-021-01101-x

[Circos: an information aesthetic for comparative genomics](https://pubmed.ncbi.nlm.nih.gov/19541911/). <b >Krzywinski M, Schein J, Birol I, Connors J, Gascoyne R, Horsman D, Jones SJ, Marra MA.</b> Genome Res. 2009 Sep;19(9):1639-45. doi:10.1101/gr.092759.109