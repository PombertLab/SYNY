## <b>Synopsis</b>

The SYNY pipeline investigates gene collinearity (synteny) between genomes by reconstructing clusters from conserved pairs of protein-coding genes identified from [DIAMOND](https://github.com/bbuchfink/diamond) homology searches. It also infers collinearity from pairwise genome alignments with [minimap2](https://github.com/lh3/minimap2).

[![DOI](https://zenodo.org/badge/491274225.svg)](https://zenodo.org/doi/10.5281/zenodo.10790180)

## <b>Table of contents</b>
* [Introduction](#Introduction)
* [Requirements](#Requirements)
  * [Downloading SYNY from GitHub](#downloading-SYNY-from-github)
  * [Installing dependencies](#installing-dependencies)
* [Using SYNY](#Using-SYNY)
  * [Why use two distinct approaches?](#Why-use-two-distinct-approaches)
  * [Command line options](#Command-line-options)
  * [Step by step examples](#Step-by-step-examples)
    * [Example 1: <i>Cryptococcus</i>](#Example-1---Cryptococcus)
      * [Heatmaps](#Heatmaps)
      * [Circos plots](#Circos-plots)
      * [Barplots (chromosome maps)](#Barplots)
      * [Dotplots](#Dotplots)
      * [PAF metrics](#PAF-metrics)
    * [Example 2: Encephalitozoonidae](#Example-2---Encephalitozoonidae)
      * [Custom Circos colors](#Custom-Circos-colors)
  * [File conversion](#File-conversion)
* [References](#References)

## <b>Introduction</b>
#### <b>What is synteny?</b>
Synteny is the measure of similarity in gene organization between two species. Through time, genomes can be reorganized by large-scale events such as recombination, but also smaller events like -- but not exhaustively -- duplication, translocation, deletion, or inversion.

#### <b>Why use synteny?</b>
Synteny inferences can be used to:
- Improve genome annotation (by helping to identify mispredicted and/or missing genes)
- Help differentiate between orthologous and paralogous genes (by leveraging positional information)
- Perform evolutionary distance analyses (<i>e.g.</i> species sharing similar/identical genome reorganization events are likely more closely related than species that do not share these events).

## <b>Requirements</b>
- [DIAMOND](https://github.com/bbuchfink/diamond)
- [minimap2](https://github.com/lh3/minimap2)
- [Perl5](https://www.perl.org/)
- [PerlIO::gzip](https://metacpan.org/pod/PerlIO::gzip)
- [Python3](https://www.python.org/)
- [matplotlib](https://matplotlib.org/)
- [seaborn](https://seaborn.pydata.org/)
- [pandas](https://pandas.pydata.org/)
- [scipy](https://scipy.org/)
- [Circos](https://circos.ca/)

#### <b>Downloading SYNY from GitHub</b>
```Bash
git clone https://github.com/PombertLab/SYNY.git
cd SYNY
export PATH=$PATH:$(pwd)
```

#### <b>Installing dependencies</b>
##### <b>Installing dependencies automatically with setup_syny.pl</b>

Dependencies can be installed automatically with `setup_syny.pl`. This script will download and install [DIAMOND](https://github.com/bbuchfink/diamond), [minimap2](https://github.com/lh3/minimap2), [Circos](https://circos.ca/) together with the required dnf/apt/zypper packages (this script has been tested on Fedora, Ubuntu, Debian, Kali, and openSUSE Tumbleweed distributions). Note that using this script will require sudo privileges to install dnf/apt/zypper packages.

```Bash
## Listing supported Linux distributions:
setup_syny.pl --list_distro

Supported Linux distributions/package managers are:

debian      apt
fedora      dnf
kali        apt
opensuse    zypper
ubuntu      apt
```

``` Bash
## Installing dependencies:

LINUX=fedora              ## Replace by Linux distribution: Ubuntu, Debian or Fedora 
CONFIG=~/.bash_profile    ## Replace by desired configuration file
DIR=~/TOOLS               ## Replace by desired installation directory

setup_syny.pl \
  --linux $LINUX \
  --config $CONFIG \
  --install $DIR

## Loading the configuration file:
source $CONFIG
```

##### <b>Installing dependencies manually</b>
##### To install PerlIO::gzip:

```Bash
## On Ubuntu:
sudo apt install libperlio-gzip-perl

## On Fedora:
sudo dnf install perl-PerlIO-gzip
```


##### To install matplolib / seaborn / pandas / scipy:
```Bash
## On Ubuntu:
sudo apt install python3-matplotlib
sudo apt install python3-seaborn
sudo apt install python3-pandas
sudo apt install python3-scipy

## On Fedora:
sudo dnf install python3-matplotlib
sudo dnf install python3-seaborn
sudo dnf install python3-pandas
sudo dnf install python3-scipy

## Or via pip (Ubuntu/Fedora):
pip install matplotlib
pip install seaborn
pip install pandas
pip install scipy
```

##### To install DIAMOND:
```Bash
version=v2.1.9     ## Replace with desired DIAMOND version
DIR=/opt/diamond   ## Replace with desired installation directory
mkdir -p $DIR

curl \
  -L https://github.com/bbuchfink/diamond/releases/download/$version/diamond-linux64.tar.gz \
  -o $DIR/diamond-linux64.tar.gz

tar -zxvf $DIR/diamond-linux64.tar.gz --directory $DIR
rm $DIR/diamond-linux64.tar.gz
export PATH=$PATH:$DIR
```

##### To install minimap2:
```Bash
## On Ubuntu
sudo apt install build-essential
sudo apt install zlib1g-dev

## On Fedora
sudo dnf group install "C Development Tools and Libraries"
sudo dnf install zlib-devel

## Ubuntu/Fedora
git clone https://github.com/lh3/minimap2
cd minimap2 && make

# To install the k8 javascript shell inside the minimap2 dir
curl \
  -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 \
  -o k8-0.2.4.tar.bz2

tar -jxf k8-0.2.4.tar.bz2
cp k8-0.2.4/k8-`uname -s` k8

export PATH=$PATH:$(pwd)      ## minimap2 install directory
export PATH=$PATH:$(pwd)/misc ## minimap2 subdir containing paftools.js
```

##### To install Circos:
```Bash
## Dependencies (Ubuntu)
sudo apt install cpanminus
sudo apt install libgd-perl

## Dependencies (Fedora)
sudo dnf install perl-App-cpanminus
sudo dnf install perl-GD

## Dependencies (Ubuntu/Fedora)
sudo cpanm \
  Clone \
  Config::General \
  Font::TTF::Font \
  List::MoreUtils \
  Math::Bezier \
  Math::Round \
  Math::VecStat \
  Params::Validate \
  Readonly \
  Regexp::Common \
  Set::IntSpan \
  Statistics::Basic \
  SVG \
  Text::Format

## Circos
wget \
  --no-check-certificate \
  https://circos.ca/distribution/circos-0.69-9.tgz

tar -zxvf circos-0.69-9.tgz
cd circos-0.69-9/bin

printf "\nexport PATH=\$PATH:$(pwd)" >> ~/.profile      ## Ubuntu
printf "\nexport PATH=\$PATH:$(pwd)" >> ~/.bash_profile ## Fedora
```

## <b>Using SYNY</b>

The SYNY pipeline can be run with [run_syny.pl](https://github.com/PombertLab/SYNY/blob/main/run_syny.pl), a master script that:
1. Extracts genome/protein sequences and annotation data from GenBank (.gbf/.gbff) flat files.
2. Performs round-robin pairwise genome alignments with [minimap2](https://github.com/lh3/minimap2).
3. Performs round-robin [DIAMOND](https://github.com/bbuchfink/diamond) BLASTP homology searches, identifies conserved protein gene pairs, and reconstructs collinear clusters from these searches.
4. Generates dotplots, barplots and [Circos](https://circos.ca/) plots highlighting collinear regions inferred from pairwise genome alignments and from shared protein cluster reconstructions.
5. Generates heatmaps summarizing the percentages of collinear bases found in pairwise genome alignments and the percentages of proteins found in collinear clusters between all genomes.

#### Why use two distinct approaches?

When working with well-annotated, closely related genomes, collinear regions inferred from pairwise genome alignments and shared protein cluster reconstructions should yield similar results.

However, when working with genomes featuring high levels of sequence divergence, pairwise genome alignments may struggle. In those instances, collinear regions inferred from protein cluster reconstructions should outperform those from genome alignments; hypervariable intergenic regions are not considered in protein cluster reconstructions and silent mutations/codon usage biases do not affect amino acid sequences. Amino acids also have a larger character state and lower back-mutation probabilities than nucleotide sequences, allowing for searches across wider evolutionary distances.

Conversely, when working with poorly annotated or unannotated genomes, collinear regions inferred from genome alignments should outperform those inferred from protein cluster reconstructions (if there is sufficient sequence similarity to perform pairwise alignments).

### Command line options

SYNY can be run from the master script as follows:<br>
```Bash
run_SYNY.pl \
  -a *.gbff \
  -o SYNY
```
Options for run_SYNY.pl are:
```
-h (--help)             Display all command line options
-t (--threads)          Number of threads to use [Default: 16]
-a (--annot)            GenBank GBF/GBFF Annotation files (GZIP files are supported)
-o (--outdir)           Output directory [Default = SYNY]
-e (--evalue)           DIAMOND BLASTP evalue cutoff [Default = 1e-10]
-g (--gaps)             Allowable number of gaps between gene pairs [Default = 0]
--asm                   Specify minimap2 max divergence preset (--asm 5, 10 or 20) [Default: off]
--resume                Resume minimap2 computations (skip completed alignments)
--no_map                Skip minimap2 pairwise genome alignments
--no_clus               Skip gene cluster reconstructions

### Circos plots
-c (--circos)           Circos plot mode: pair (pairwise), cat (concatenated), all (cat + pair) [Default: pair]
--orientation           Karyotype orientation: normal, inverted or both [Default: normal]
--circos_prefix         Desired Circos plot prefix for concatenated plots [Default: circos]
-r (--ref)              Genome to use as reference for concatenated plots (defaults to first one alphabetically if none provided)
-u (--unit)             Size unit (Kb or Mb) [Default: Mb]
--winsize               Sliding windows size (nucleotide biases) [Default: 10000]
--stepsize              Sliding windows step (nucleotide biases) [Default: 5000]
--labels                Contig label type: numbers or names [Default: numbers]
--label_size            Contig label size [Default: 36]
--label_font            Contig label font [Default: bold] - https://circos.ca/documentation/tutorials/ideograms/labels/
--custom_file           Load custom colors from file
--list_preset           List available custom color presets
--custom_preset         Use a custom color preset, e.g.: --custom_preset chloropicon
--max_ticks             Set max number of ticks [Default: 5000]
--max_ideograms         Set max number of ideograms [Default: 200]
--max_links             Set max number of links [Default: 25000]
--max_points_per_track  Set max number of points per track [Default: 75000]
--clusters              Color by cluster instead of contig/chromosome [Default: off]
--no_circos             Turn off Circos plots

### Barplots
-bh (--bheight)         Barplot figure height in inches [Default: 10.8]
-bw (--bwidth)          Barplot figure width in inches [Default: 19.2]
--palette               Barplot color palette [Default: Spectral]
--monobar               Use a monochrome barplot color instead: e.g. --monobar blue
--no_barplot            Turn off barplots

### Dotplots
-dh (--dheight)         Dotplot figure height in inches [Default: 10.8]
-dw (--dwidth)          Dotplot figure width in inches [Default: 19.2]
-m (--multi)            Axes units multiplier (for dotplots) [Default: 1e5]
--color                 Dotplot color [Default: blue]
--dotpalette            Use a color palette instead: e.g. --dotpalette inferno
--noticks               Turn off ticks on x and y axes
--wdis                  Horizontal distance (width) between subplots [Default: 0.05]
--hdis                  Vertical distance (height) between subplots [Default: 0.1]
--no_dotplot            Turn off dotplots

### Heatmaps
-hh (--hheight)         Heatmap figure height in inches [Default: 10]
-hw (--hwidth)          Heatmap figure width in inches [Default: 10]
--hmpalette             Heatmap color palette [Default: winter_r]
--no_heatmap            Turn off heatmaps
```
The output directory will be structured as follows: 
```Bash
ls -lah SYNY/

drwxr-xr-x  7 jpombert jpombert 4096 Apr 24 09:54 .
drwx------ 20 jpombert jpombert 4096 Apr 24 09:54 ..
drwxr-xr-x  6 jpombert jpombert 4096 Apr 24 09:54 ALIGNMENTS
drwxr-xr-x  6 jpombert jpombert 4096 Apr 24 09:55 CLUSTERS
drwxr-xr-x  2 jpombert jpombert 4096 Apr 24 09:54 LISTS
drwxr-xr-x  7 jpombert jpombert 4096 Apr 24 09:55 PLOTS
drwxr-xr-x  4 jpombert jpombert 4096 Apr 24 09:54 SEQUENCES
-rw-r--r--  1 jpombert jpombert  526 Apr 24 09:55 error.log
-rw-r--r--  1 jpombert jpombert 1574 Apr 24 09:57 syny.log
```

Genome/protein sequences and gene lists extracted from GenBank GBFF files are located in the `SEQUENCES/` and `LISTS/` directories, respectively. Data files from minimap2 pairwise genome alignments and from gene cluster inferences are located in the `ALIGNMENTS/` and `CLUSTERS/` directories, respectively. Plots are located in `PLOTS/`.

The contents of the subdirectories are:
- ALIGNMENTS:
	- MAF, PAF and ALN: pairwise genome alignments in the corresponding formats
	- METRICS:
      - Alignment length vs. similarity scatter plots (in PNG format)
      - Alignment metrics summaries (in plain TXT format)
- CLUSTERS:
  - DIAMOND:
    - DB
      - BLASTP databases
    - Round-robin BLASTP results (.diamond.6)
  - HOMOLOGS:
    - Lists of proteins that are conserved within the samples (.conserved)
    - Homology and conservation summaries (.conserved_summary)
    - Lists of unique proteins and their locations (.unique)
    - Tab-delimited protein/description lists of proteins and their top homologs in other species (.shared.tsv)
    - Tab-delimited protein/description lists of proteins uniques to each species (.uniques.tsv)
    - Tab-delimited protein/description lists of all proteins, shared or not (.all.tsv)
  - SYNTENY:
    - Cluster summary (clusters_summary.tsv)
    - Tab-delimited cluster summary table (clusters_summary_table.tsv)
    - Subdirectory per specified gap allowance (gap_#) containing:
      - CLUSTERS:
        - Round-robin reconstructed syntenic clusters for each species
      - PAIRS:
        - Round-robin identified gene pairs for each species
- LISTS:
	- Lists of protein coding genes with location details (.list)
- PLOTS:
  - BARPLOTS
	  - Barplots (in PNG/SVG format) from minimap2 PAF alignments (.mmap.)
	  - Barplots (in PNG/SVG format) from protein clusters found with SYNY (.gap_0., .gap_1., ...)
  - CIRCOS:
	  - Circos plots (in PNG/SVG format) from minimap2 PAF alignments (.mmap.)
	  - Circos plots (in PNG/SVG format) from protein clusters found with SYNY (.gap_0., .gap_1., ...)
  - CIRCOS_DATA:
	  - Configuration files for Circos plots
  - DOTPLOTS:
  	- Dotplots (in PNG/SVG format) from minimap2 PAF alignments (.mmap.)
  	- Dotplots (in PNG/SVG format) from protein clusters found with SYNY (.gap_0., .gap_1., ...)
  - HEATMAPS:
    - Heatmaps (in PNG/SVG format) summarizing the percentages of collinear bases in pairwise alignments (.mmap.)
    - Heatmaps (in PNG/SVG format) summarizing the percentages of proteins found in clusters (.gap_0., .gap_1., ...)
- SEQUENCES:
  - GENOMES:
  	- FASTA files containing the sequences of the investigated genomes
  - PROTEINS:
    - Protein sequences for each species (.faa) in FASTA format

### <b>Step by step examples</b>

### Example 1 - <i>Cryptococcus</i>
Below is a quick example describing how to compare a few select <i>Cryptococcus</i> genomes (<i>C. neoformans</i> var. <i>neoformans</i> strain [JEC21](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA13856/), <i>C. gattii</i> strain [WM276](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA13692/), <i>C. gattii</i> VGV strain [MF34](hhttps://www.ncbi.nlm.nih.gov/bioproject/PRJNA487802/), <i>C. decagattii</i> strain [7685027](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA660466/), <i>C. deuterogattii</i> strain [R265](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA395628/)) using annotation data available from the public NCBI [Genome](https://www.ncbi.nlm.nih.gov/genome/database) database.

##### Downloading annotation data from GenBank (NCBI):

Downloading <i>Cryptococcus</i> example data automatically from NCBI using the provided [Cryptococcus.sh](https://github.com/PombertLab/SYNY/blob/main/Examples/Cryptococcus.sh) shell script from the `Examples/` subdirectory:
```bash
DATA=~/DATA ## Replace ~/DATA by desired annotation data directory
Cryptococcus.sh $DATA
```

Or downloading the <i>Cryptococcus</i> data manually instead:
```bash
DATA=~/DATA ## Replace ~/DATA by desired annotation data directory
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
```

##### Running SYNY:

SYNY can be run from the command line with the `run_syny.pl` master script. In the command line below, no gap is allowed during gene cluster inferences, and Circos plots are produced in pairwise mode. The latter can be changed with the `--circos` command line switch (possible values are: `pair`, `concatenated`, `both`). Note that when comparing several genomes, concatenated plots can quickly become too dense for legibility. Producing concatenated plots can also significantly increase computation time.

```Bash
##### Runtime (Intel i5-12500H mobile CPU)
# Total: 16 min; Circos: 11 min

DATA=~/DATA                   ## Replace by annotation data directory
SYNY=~/SYNY_RESULTS/CRYPT_ALL ## Replace by desired output directory

run_syny.pl \
  --threads 16 \
  --annot $DATA/*.gbff.gz \
  --outdir $SYNY
```

In the command line below, SYNY is run on a subset of two genomes (`JEC21` and `WM276`), this time allowing for different maximum gap thresholds between gene pairs.

```Bash
##### Runtime (Intel i5-12500H mobile CPU)
# Total: 2 min; Circos: 1 min)

DATA=~/DATA                   ## Replace by annotation data directory
SYNY=~/SYNY_RESULTS/CRYPT_SUB ## Replace by desired output directory

run_syny.pl \
  -t 16 \
  -a $DATA/{JEC21,WM276}.gbff.gz \
  -o $SYNY \
  -g 0 1 5 \
  -e 1e-10
```

##### Example of clusters identified with SYNY
```Bash
head -n 32 $SYNY/CLUSTERS/SYNTENY/clusters_summary.tsv

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
head -n 22 $SYNY/CLUSTERS/SYNTENY/gap_0/CLUSTERS/JEC21_vs_WM276.clusters

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

Overall cluster metrics are summarized in `CLUSTERS/SYNTENY/clusters_summary_table.tsv`.

```Bash
head -n 7 $SYNY/CLUSTERS/SYNTENY/clusters_summary_table.tsv

### Query       Total # proteins        Allowed Gaps    Total # proteins in clusters    % of proteins in clusters       # of clusters   Longest Shortest        Average  Median  N50     N75     N90
JEC21_vs_WM276  6863    0       5762    83.96   779     64      2       7       5       10      6       4
JEC21_vs_WM276  6863    1       5922    86.29   228     204     2       26      15      53      26      13
JEC21_vs_WM276  6863    5       5957    86.80   53      434     2       112     58      237     163     66
WM276_vs_JEC21  6565    0       5760    87.74   778     64      2       7       5       10      6       4
WM276_vs_JEC21  6565    1       5923    90.22   229     204     2       26      15      53      26      13
WM276_vs_JEC21  6565    5       5957    90.74   51      477     2       117     56      260     164     66
```

#### Heatmaps

To facilitate comparisons when working with large datasets, heatmaps displaying the percentages of collinear bases in pairwise genome alignments (`.mmap.`) and the percentages of protein coding-genes found in collinear clusters between each pair of genomes (e.g. `.gap_0.`) are generated with matplotlib.

In the above example, small heatmaps with 25 datapoints (5 * 5 genomes) will be generated for the minimap2 pairwise alignments and for the protein clusters found for each gap value investigated. In these heatmaps, percentages between pairs of genomes will vary based on the total number of bases and proteins found in the query used: <i>i.e.</i> `(collinear bases / total bases) * 100` and `(proteins in clusters / total proteins) * 100`.

##### Example of a heatmap showing percentages of collinear protein-coding genes between investigated genomes:
<p align="left">
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/proteins_in_clusters.gap_0.heatmap.png">
</p>

Heatmap dimensions (default: 10 x 10) can be modified with the `--hheight` and `--hwidth` command line switches. The color palette (default: winter_r) can be modified with the `--hmpalette` command line switch (see this [URL](https://www.practicalpythonfordatascience.com/ap_seaborn_palette) for a list of available seaborn color palettes).

#### Circos plots

Unless the `--no_circos` command line switch is invoked, [Circos](https://circos.ca/) (pairwise and/or concatenated) will be generated from the protein clusters identified with SYNY (e.g. `.gap_0.`) and/or from the genome alignments computed with minimap2 (`.mmap.`).

In the pairwise plots (`--circos pair`), genomes are plotted in pairs (query <i>vs.</i> subject) using the query as the reference. In the concatenated plots (`--circos cat`), all genomes are plotted together in a single figure, using the reference genome specified with the `--ref` command line switch. If omitted, the first genome encountered alphabetically will be used as the default reference. Both concatenated and pairwise plots can be generated with the `--circos all` command line switch.

Karyotypes can be plotted with Circos in normal and/or inverted orientation(s). In the inverted plots, the contigs/chromosomes from the genome(s) being compared to the reference are plotted in reverse, from last to first. By default, the Circos plots are plotted in normal orientation. This behavior can be modified with the `--orientation` command line switch (possible values are `normal`, `inverted` and `both`).

##### Example of an image generated with Circos and SYNY from shared proteins clusters:
<p align="left">
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/WM276_vs_JEC21.gap_0.normal.png">
</p>

In this figure, nucleotides biases are plotted in the concentric rings (from outer to inner rings):
- AT and GC nucleotide biases (grey and red lines)
- GT and AC nucleotide biases (blue and green lines)
- GA and CT nucleotide biases (purple and yellow lines)

Syntenic blocks identified by SYNY are indicated by ribbons. These ribbons are color-coded based on the chromosomes/contigs present in the reference genome used. If desired, color coding can be set by cluster instead with `--clusters`. This option is useful when working with prokaryotes featuring a single chromosome, so that the ribbons are not all of the same color.

By default, contigs will be labelled by numbers in the Circos plots. If desired, contigs can instead be labelled by their names with the `--labels names` command line option. Label sizes (default: 36) and fonts can be further adjusted with the `--label_size` and `--label_font` command line options. Possible fonts are: `light`, 
`normal`, `default`, `semibold`, `bold`, `italic`, `bolditalic`, `italicbold` (see this Circos [tutorial](https://circos.ca/documentation/tutorials/ideograms/labels/) for details).

```Bash
##### Runtime (Intel i5-12500H mobile CPU)
# Total: 2 min; Circos: 1 min

DATA=~/DATA                         ## Replace by annotation data directory
SYNY=~/SYNY_RESULTS/CRYPT_SUB_NAMES ## Replace by desired output directory

run_syny.pl \
  -a $DATA/{JEC21,WM276}.gbff.gz \
  -g 0 1 5 \
  -e 1e-10 \
  -r JEC21 \
  -o $SYNY \
  --circos all \
  --circos_prefix WM276_vs_JEC21 \
  --labels names \
  --label_size 20 \
  --label_font semibold
```

##### Example of the same Circos plot, this time labelled by contig names intead:
<p align="left">
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/WM276_vs_JEC21.gap_0.normal.names.png">
</p>

All data and configuration files for the Circos plots are in the `PLOTS/CIRCOS_DATA/` subdirectory. If desired, chromosome labels can be further adjusted manually by editing the `LABEL` column in the corresponding Circos [karyotype](https://circos.ca/documentation/tutorials/ideograms/karyotypes/) file (`concatenated.normal.genotype` / `concatenated.inverted.genotype` from `CIRCOS/concatenated/`).

```Bash
## To see the first 5 lines from a karyotype file:
head -n 5 concatenated.normal.genotype

#chr - ID LABEL START END COLOR
chr - NC_014938 1 0 1984822 black
chr - NC_014939 2 0 2187694 black
chr - NC_014940 3 0 1961511 black
chr - NC_014941 4 0 2233617 black
```

Once edited, Circos plots can be regenerated by running Circos on the corresponding configuration files, <i>e.g.</i>:

```Bash
SYNY=~/SYNY_RESULTS/CRYPT_SUB_NAMES  ## Replace by SYNY results directory
OUTDIR=~/PLOTS                       ## Replace by desired output directory

circos \
  -conf $SYNY/PLOTS/CIRCOS_DATA/concatenated/concatenated.gap_0.normal.conf \
  -outputdir $OUTDIR \
  -outputfile circos.gap_0.normal.png
```

##### Example of an image generated with Circos and SYNY from shared proteins clusters (using the inverted karyotype):
<p align="left">
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/WM276_vs_JEC21.gap_0.inverted.png">
</p>

In the inverted karyotype image, the order of the karyotype(s) to be compared to the reference one is reversed. This option can be useful when comparing genomes whose chromosomes have been assigned similar numbers based on various inference methods (this does not appear to be the case in the above example). In such instances, inverting the karyotypes can help improve figure legibility.

##### Example of an image generated with Circos and SYNY from pairwise genome alignments:
<p align="left">
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/WM276_vs_JEC21.mmap.normal.png">
</p>

In pairwise genome alignments, repetitive regions (such as telomeres/subtelomeres) can produce more than one alignment for a given locus. In the above figure, a bit of extra noise is added to the figure (as thin criss-crossing lines) due to these repetitive segments. As a rule of thumb, repetitive segments are easier to spot in dotplot-like figures (see [Dotplots](https://github.com/PombertLab/SYNY?tab=readme-ov-file#Dotplots) section below).

#### Barplots
Chromosome maps (aka barplots) highlighting collinear genome segments will be generated with [paf_to_barplot.py](https://github.com/PombertLab/SYNY/blob/main/Plots/paf_to_barplot.py) and [matplotlib](https://matplotlib.org/)  from [minimap2](https://github.com/lh3/minimap2) pairwise genome alignments (`.mmap.`) and from protein clusters identified with SYNY (e.g. `.gap_0.`).

##### Example of a barplot generated from minimap2 PAF files using defaults settings:
<p align="left">
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/WM276_vs_JEC21.mmap.barplot.19.2x10.8.Spectral.png">
</p>

In these plots, collinear regions found between the compared genomes are highlighted by colored rectangles. By default, these rectangles are color-coded based on the contigs/chromosomes of the query. The above barplot image was generated using the Spectral color palette from [seaborn](https://seaborn.pydata.org/tutorial/color_palettes.html) (set as default in SYNY). This palette can be replaced using the `--palette` command line switch; <i>e.g.</i> `--palette husl` (see this [URL](https://www.practicalpythonfordatascience.com/ap_seaborn_palette) for a detailed list of available palettes).

If desired, the barplots can instead be generated using a single monochromatic color with the `--monobar` command line switch; <i>e.g.</i> `--monobar red`. If this option is selected, the color-coding legend will be omitted from the plot.

By default, the barplots are formatted for a widescreen (landscape) output (width/height ratio: 19.2/10.8). This ratio can be adjusted with the `--bheight` and `--bwidth` command line switches.

#### Dotplots

Unless the `--no_dotplot` command line switch is invoked, dotplots will be generated with [paf_to_dotplot.py](https://github.com/PombertLab/SYNY/blob/main/Plots/paf_to_dotplot.py) and [matplotlib](https://matplotlib.org/) from [minimap2](https://github.com/lh3/minimap2) pairwise genome alignments (`.mmap.`) and from protein clusters identified with SYNY, (e.g. `.gap_0.`).

##### Example of a dotplot generated from minimap2 PAF files using defaults settings:

<p align="left">
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/WM276_vs_JEC21.mmap.1e5.19.2x10.8.blue.png">
</p>

In these dotplots, each chromosome/contig from the query is plotted as a column (x-axis) against each chromosome/contig from the subject (y-axis). In the above example, a total of 196 subplots (14 x 14 chromosomes) are plotted using matplotlib's [subplot](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.subplots.html) function. In these plots, matches are plotted using the selected color (defaut: blue), with repeated loci indicated by the presence of matches across two or more contigs/chromosomes.

Because large numbers quickly overlap in the small subplots, to improve legibility, numbers in the x and y axes are scaled down using a desired scale with the `--multi` command line switch (default: 1e5). For example:
```bash
1e3: n x 1000 bp
1e4: n x 10000 bp
1e5: n x 100000 bp
1e6: n x 1000000 bp
```

Alternatively, ticks/numbers in the x and y axes can be turned off with the `--noticks` command line option.

By default, the dotplots are formatted for a widescreen (landscape) output (width/height ratio: 19.2/10.8). This ratio can be adjusted with the `--dheight` and `--dwidth` command line switches.

The default monochromatic color (blue) can be changed with the `--color` option (<i>e.g.</i> `--color red`). If desired, dotplots can instead be color-coded based on the query contigs/chromosomes with the `--dotpalette` option; <i>e.g.</i> `--dotpalette husl` (see this [URL](https://www.practicalpythonfordatascience.com/ap_seaborn_palette) for a detailed list of available seaborn palettes).

##### Example of a dotplot generated from minimap2 PAF files using the husl color palette:

```Bash
##### Runtime (Intel i5-12500H mobile CPU);
# Total: < 1 min; Circos: N/A

DATA=~/DATA                     ## Replace by annotation data directory
SYNY=~/SYNY_RESULTS/CRYPT_HUSL  ## Replace by desired output directory

## Running SYNY with the husl color palette (for dotplots)
## and the --no_circos option to skip Circos plotting
run_syny.pl \
  -a $DATA/{JEC21,WM276}.gbff.gz \
  -g 0 1 5 \
  -e 1e-10 \
  -o $SYNY \
  --dotpalette husl \
  --no_circos
```

<p align="left">
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/WM276_vs_JEC21.mmap.1e5.19.2x10.8.husl.png">
</p>

#### PAF metrics

When computing pairwise genome alignments with [minimap2](https://github.com/lh3/minimap2), PAF metrics will also be calculated independently to help assess the outcome of these alignments. These metrics are available in the `ALIGNMENTS/METRICS/` subdirectory, with alignment length <i>vs.</i> sequence similarity (%) scatter plots available in PNG format.

##### Example of a scatter plot summarizing metrics from minimap2 PAF files:

<p align="left">
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/WM276_vs_JEC21.metrics.png">
</p>

In the above plot, the average sequence identity for each alignment is calculated from the [PAF](https://github.com/lh3/miniasm/blob/master/PAF.md) files as follows:

$(\\#\ of\ residue\ matches / Alignment\ block\ length) * 100$

Likewise, the total average sequence identity percentage listed in the PAF metrics insert is calculated by dividing the sum of all residue matches by the sum of all aligned block lengths.

As a rule of thumb, pairwise alignments featuring lower sequence identity percentages will produce fewer and/or more fragmented collinear segments.

### Example 2 - Encephalitozoonidae
Below is a quick example describing how to compare a few select genomes from the Encephalitozoonidae (<i>Encephalitozoon/Ordospora</i> species <i>E. intestinalis</i> strain [ATCC 50506](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA594722/), <i>E. hellem</i> strain [ATCC 50604](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA594722/), <i>E. cuniculi</i> strain [ATCC 50602](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA705735/), <i>O. colligata</i> strain [OC4](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA210314/) and <i>O. pajunii</i> strain [FI-F-10](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA630072/)) using annotation data from NCBI.

The Encephalitozoonidae genomes in this example are highly collinear yet display a high level of sequence divergence. As such, they constitute a good case scenario to assess the benefits of using gene clusters to complement inferences based on pairwise genome alignments.

##### Downloading annotation data from GenBank (NCBI):

Downloading genome data automatically from NCBI using the [Encephalitozoonidae.sh](https://github.com/PombertLab/SYNY/blob/main/Examples/Encephalitozoonidae.sh) shell script from the `Examples/` subdirectory:

```bash
DATA=~/ENCE ## Replace ~/ENCE by desired annotation data directory
Encephalitozoonidae.sh $DATA
```
Or downloading the Encephalitozoonidae data manually instead:

```bash
DATA=~/ENCE ## Replace ~/ENCE by desired annotation data directory
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
```

##### Running SYNY:
To reduce runtime, we will skip the Circos plots with the `--no_circos` option.

```Bash
##### Runtime (Intel i5-12500H mobile CPU);
# Total: 1 min; Circos: N/A

DATA=~/ENCE                ## Replace by annotation data directory
SYNY=~/SYNY_RESULTS/ENCE   ## Replace by desired output directory

run_syny.pl \
  -a $DATA/*.gbff.gz \
  -g 0 \
  -e 1e-10 \
  -o $SYNY \
  --no_circos
```

Here, while the pairwise alignments and gene cluster collinear inferences produce similar results within species (not shown), they differ substantially between species. For example, below are two barplots between <i>Encephalitozoon hellem</i> and <i>Ordospora colligata</i>. While the results are congruent, gene clusters inferences `.gap_0.` perform better than those based on pairwise alignments (`.mmap.`), the latter of which struggle due to the high levels of sequence divergence involved. These differences will also be reflected in the dotplots, heatmaps and Circos plots generated (not shown).

Barplots (gene clusters): 
<p align="left">
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/Eh50604_vs_OcOC4.gap_0.barplot.19.2x10.8.Spectral.png">
</p>

Barplots (genome alignments): 
<p align="left">
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/Eh50604_vs_OcOC4.mmap.barplot.19.2x10.8.Spectral.png">
</p>

##### About concatenated Circos plots:
Because comparing more than two genomes can be useful, concatenated Circos plots can be generated by SYNY using the `--circos all` or the` --circos cat` command line switches. Below is a quick example of how to generate such plots using the genomes from the three <i>Encephalitozoon</i> species downloaded previously.

```Bash
##### Runtime (Intel i5-12500H mobile CPU);
# Total: < 1 min; Circos: 17 sec

DATA=~/ENCE                    ## Replace by annotation data directory
SYNY=~/SYNY_RESULTS/SYNY_3_spp ## Replace by desired output directory

run_syny.pl \
  -a $DATA/E*.gbff.gz \
  -g 0 \
  -e 1e-10 \
  -r Ei50506 \
  -o $SYNY \
  --circos cat \
  --circos_prefix encephalitozoon
```


##### Example of a concatenated Circos plot comparing a total of 3 genomes (using defaults settings):
<p align="left">
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/encephalitozoon.syny.normal.png">
</p>

In the above [Circos](https://circos.ca/) plot, links between the reference and the queried genomes are color-coded based on the reference genotype; links between the queried genomes are shown in light gray. In these plots, each chromosome/contig is plotted as a distinct [ideogram](https://circos.ca/documentation/tutorials/ideograms/ideograms/).

Note that while SYNY can plot multiple genomes together with `--circos cat` / `--circos all`, adding too much data can quickly clutter the concatenated Circos plots (drawing too many ideograms can make the figure illegible). If the concatenated plots are too dense, drawing pairwise plots with `--circos pair` should help with readability. However, when comparing genomes featuring many chromosomes/contigs, cluttering can occur even in pairwise plots.

#### Custom Circos colors
Custom colors for [Circos](https://circos.ca/) plots can be loaded directly from tab-delimited text files containing color names and their associated RGB values (see [custom_color_1.txt](https://github.com/PombertLab/SYNY/blob/main/Examples/custom_color_1.txt) for an example). A few custom presets are also available to use in SYNY.

##### Running SYNY with custom colors loaded from a tab-delimited file:
```Bash
##### Runtime (Intel i5-12500H mobile CPU);
# Total: 1 min; Circos: 37 sec

DATA=~/ENCE                  ## Replace by annotation data directory
SYNY=~/SYNY_RESULTS/ENCE_CC  ## Replace by desired output directory
COLORS=~/custom_color_2.txt  ## Replace by desired custom color file

run_syny.pl \
  -a $DATA/E*.gbff.gz \
  -g 0 1 5 \
  -e 1e-10 \
  -r Ei50506 \
  -o $SYNY \
  --custom_file $COLORS \
  --circos all \
  --circos_prefix encephalitozoon_cc
```

##### Example of an image generated with Circos and SYNY comparing a total of 3 genomes (using a custom color set):
<p align="left">
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/encephalitozoon_cc.syny.normal.png">
</p>

##### Listing custom color presets available in SYNY:
```
run_syny.pl --list_preset

Available color presets:
blues   13 colors
chloropicon     20 colors
encephalitozoon 11 colors
```

##### Running SYNY with a custom color preset:
```Bash
##### Runtime (Intel i5-12500H mobile CPU);
# Total: 1 min; Circos: 37 sec

DATA=~/ENCE                        ## Replace by annotation data directory
SYNY=~/SYNY_RESULTS/ENCE_CC_PRESET ## Replace by desired output directory

run_syny.pl \
  -a $DATA/E*.gbff.gz \
  -g 0 1 5 \
  -e 1e-10 \
  -r Ei50506 \
  -o $SYNY \
  --circos all \
  --circos_prefix encephalitozoon_blues \
  --custom_preset blues
```

##### Example of an image generated with the blues custom color preset:
<p align="left">
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/encephalitozoon_blues.syny.normal.png">
</p>

### <b>File conversion</b>
SYNY uses [GenBank flat files](https://www.ncbi.nlm.nih.gov/genbank/samplerecord/) (.gbff) containing DNA sequences and annotated features as input. However, when performing genome alignment-based comparisons, the latter annotated features are not required. Because users may want to compare genomes before their annotation becomes available (this process often comes in the late stages of genome sequencing projects), a simple FASTA to GBFF converter [fasta_to_gbff.pl](https://github.com/PombertLab/SYNY/blob/main/Utils/fasta_to_gbff.pl) producing feature-less .gbff files is available in the `Utils/` subdirectory.

To convert FASTA file(s) to feature-less GBFF files compressed in gzip format:
```Bash
fasta_to_gbff.pl \
  --fasta *.fasta.gz \
  --outdir GBFF \
  --gzip
```

Options for `fasta_to_gbff.pl` are:
```
-f (--fasta)    FASTA file(s) to convert (gziped files are supported)
-o (--outdir)   Output directory [Default: GBFF]
-g (--gzip)     Compress the GBFF output files
-v (--verbose)  Add verbosity
```
When compressing the GBFF output files, `fasta_to_gbff.pl` will use `pigz` if available, otherwise it will default to `gzip`.

## <b>References</b>
[Sensitive protein alignments at tree-of-life scale using DIAMOND](https://www.nature.com/articles/s41592-021-01101-x). Buchfink B, Reuter K, Drost HG. <b>Nature Methods.</b> 18, 366–368 (2021). doi: 10.1038/s41592-021-01101-x

[Minimap2: pairwise alignment for nucleotide sequences](https://pubmed.ncbi.nlm.nih.gov/29750242/). Li H. <b>Bioinformatics.</b> 2018 Sep 15;34(18):3094-3100. doi: 10.1093/bioinformatics/bty191.

[Circos: an information aesthetic for comparative genomics](https://pubmed.ncbi.nlm.nih.gov/19541911/). Krzywinski M, Schein J, Birol I, Connors J, Gascoyne R, Horsman D, Jones SJ, Marra MA. <b>Genome Res.</b> 2009 Sep;19(9):1639-45. doi: 10.1101/gr.092759.109.

[Matplotlib: A 2D graphics environment](https://doi.org/10.1109/MCSE.2007.55). Hunter, JD. <b>Computing in Science & Engineering.</b> 2007 9(3):90-95. doi: 10.1109/MCSE.2007.55.

[Seaborn: statistical data visualization](https://doi.org/10.21105/joss.03021). Waskom, ML. <b>Journal of Open Source Software.</b> 2021 6(60): 3021. doi: 10.21105/joss.03021

[Data structures for statistical computing in Python](https://doi.org/10.25080/Majora-92bf1922-00a). McKinney W. <b>Proceedings of the 9th Python in Science Conference.</b> 2010:56-61. doi: 10.25080/Majora-92bf1922-00a.

[SciPy 1.0: fundamental algorithms for scientific computing in Python.](https://pubmed.ncbi.nlm.nih.gov/32015543/). Virtanen P, Gommers R, Oliphant TE, Haberland M, Reddy T, Cournapeau D, Burovski E, Peterson P, Weckesser W, Bright J, van der Walt SJ, Brett M, Wilson J, Millman KJ, Mayorov N, Nelson ARJ, Jones E, Kern R, Larson E, Carey CJ, Polat İ, Feng Y, Moore EW, VanderPlas J, Laxalde D, Perktold J, Cimrman R, Henriksen I, Quintero EA, Harris CR, Archibald AM, Ribeiro AH, Pedregosa F, van Mulbregt P; SciPy 1.0 Contributors. <b>Nat Methods.</b> 2020 Mar;17(3):261-272. doi: 10.1038/s41592-019-0686-2. Epub 2020 Feb 3. PMID: 32015543; PMCID: PMC7056644.

[The genome of the basidiomycetous yeast and human pathogen <i>Cryptococcus neoformans</i>](https://pubmed.ncbi.nlm.nih.gov/15653466/). Loftus BJ, Fung E, Roncaglia P, Rowley D, Amedeo P, Bruno D, Vamathevan J, Miranda M, Anderson IJ, Fraser JA, Allen JE, Bosdet IE, Brent MR, Chiu R, Doering TL, Donlin MJ, D'Souza CA, Fox DS, Grinberg V, Fu J, Fukushima M, Haas BJ, Huang JC, Janbon G, Jones SJ, Koo HL, Krzywinski MI, Kwon-Chung JK, Lengeler KB, Maiti R, Marra MA, Marra RE, Mathewson CA, Mitchell TG, Pertea M, Riggs FR, Salzberg SL, Schein JE, Shvartsbeyn A, Shin H, Shumway M, Specht CA, Suh BB, Tenney A, Utterback TR, Wickes BL, Wortman JR, Wye NH, Kronstad JW, Lodge JK, Heitman J, Davis RW, Fraser CM, Hyman RW. <b>Science.</b> 2005 Feb 25;307(5713):1321-4. doi: 10.1126/science.1103773. Epub 2005 Jan 13. PMID: 15653466; PMCID: PMC3520129.

[Genome variation in <i>Cryptococcus gattii</i>, an emerging pathogen of immunocompetent hosts](https://pubmed.ncbi.nlm.nih.gov/21304167/). D'Souza CA, Kronstad JW, Taylor G, Warren R, Yuen M, Hu G, Jung WH, Sham A, Kidd SE, Tangen K, Lee N, Zeilmaker T, Sawkins J, McVicker G, Shah S, Gnerre S, Griggs A, Zeng Q, Bartlett K, Li W, Wang X, Heitman J, Stajich JE, Fraser JA, Meyer W, Carter D, Schein J, Krzywinski M, Kwon-Chung KJ, Varma A, Wang J, Brunham R, Fyfe M, Ouellette BF, Siddiqui A, Marra M, Jones S, Holt R, Birren BW, Galagan JE, Cuomo CA. <b>mBio</b>. 2011 Feb 8;2(1):e00342-10. doi: 10.1128/mBio.00342-10. PMID: 21304167; PMCID: PMC3037005.

[A New Lineage of Cryptococcus gattii (VGV) Discovered in the Central Zambezian Miombo Woodlands](https://pubmed.ncbi.nlm.nih.gov/31719178/). Farrer RA, Chang M, Davis MJ, van Dorp L, Yang DH, Shea T, Sewell TR, Meyer W, Balloux F, Edwards HM, Chanda D, Kwenda G, Vanhove M, Chang YC, Cuomo CA, Fisher MC, Kwon-Chung KJ. <b>mBio</b>. 2019 Nov 12;10(6):e02306-19. doi: 10.1128/mBio.02306-19. PMID: 31719178; PMCID: PMC6851281.

[Telomere-to-Telomere genome assemblies of human-infecting <i>Encephalitozoon</i> species](https://pubmed.ncbi.nlm.nih.gov/37142951/). Mascarenhas Dos Santos AC, Julian AT, Liang P, Juárez O, Pombert JF. <b>BMC Genomics</b>. 2023 May 4;24(1):237. doi: 10.1186/s12864-023-09331-3. PMID: 37142951; PMCID: PMC10158259.

[A new microsporidian parasite, Ordospora pajunii sp. nov (Ordosporidae), of Daphnia longispina highlights the value of genomic data for delineating species boundaries](https://pubmed.ncbi.nlm.nih.gov/35279911/). de Albuquerque NRM, Haag KL, Fields PD, Cabalzar A, Ben-Ami F, Pombert JF, Ebert D. <b>J Eukaryot Microbiol</b>. 2022 May;69(3):e12902. doi: 10.1111/jeu.12902. Epub 2022 Mar 28. PMID: 35279911.

[The Ordospora colligata genome: Evolution of extreme reduction in microsporidia and host-to-parasite horizontal gene transfer](https://pubmed.ncbi.nlm.nih.gov/25587016/). Pombert JF, Haag KL, Beidas S, Ebert D, Keeling PJ. <b>mBio</b>. 2015 Jan 13;6(1):e02400-14. doi: 10.1128/mBio.02400-14. PMID: 25587016; PMCID: PMC4313915.