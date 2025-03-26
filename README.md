<p align="left"><img src="https://github.com/PombertLab/SYNY/blob/main/Images/SYNY.logo.png" alt="SYNY - Genome collineariy inferences" width="800"></p>

## <b>Synopsis</b>

The SYNY pipeline investigates gene collinearity (synteny) between genomes by reconstructing clusters from conserved pairs of protein-coding genes identified from [DIAMOND](https://github.com/bbuchfink/diamond) homology searches. It also infers collinearity from pairwise genome alignments with [minimap2](https://github.com/lh3/minimap2) or [MashMap3](https://github.com/marbl/MashMap).

[![DOI](https://zenodo.org/badge/491274225.svg)](https://zenodo.org/doi/10.5281/zenodo.10790180)

<details open>
  <summary><b><i>Show/hide: TOC</i></b></summary>

## <b>Table of contents</b>
* [Introduction](#Introduction)
* [Requirements](#Requirements)
  * [Downloading SYNY from GitHub](#downloading-SYNY-from-github)
  * [Installing SYNY](#installing-SYNY)
    * [With Conda](#With-Conda)
    * [With Perl install script](#With-Perl-install-script)
    * [Manual installation](#Manual-installation)
  * [About memory usage](#about-memory-usage)
* [Using SYNY](#Using-SYNY)
  * [Why use two distinct approaches?](#Why-use-two-distinct-approaches)
  * [Command line options](#Command-line-options)
  * [Step by step examples](#Step-by-step-examples)
    * [Example 1: <i>Cryptococcus</i>](#Example-1---Cryptococcus)
      * [Heatmaps](#Heatmaps)
      * [Circos plots](#Circos-plots)
      * [Barplots (chromosome maps)](#Barplots)
      * [Linemaps (linear maps)](#Linemaps)
      * [Dotplots](#Dotplots)
      * [PAF metrics](#PAF-metrics)
    * [Example 2: Encephalitozoonidae](#Example-2---Encephalitozoonidae)
      * [Custom Circos colors](#Custom-Circos-colors)
    * [Example 3: Subsets](#Example-3---Subsets)
  * [File conversion](#File-conversion)
* [Funding and acknowledgments](#Funding-and-acknowledgments)
* [How to cite](#how-to-cite)
* [References](#References)
</details>

<details open>
  <summary><b><i>Show/hide: Introduction</i></b></summary>

## <b>Introduction</b>
#### <b>What is synteny?</b>
Synteny is the measure of similarity in gene organization between two species. Through time, genomes can be reorganized by large-scale events such as recombination, but also smaller events like -- but not exhaustively -- duplication, translocation, deletion, or inversion.

#### <b>Why use synteny?</b>
Synteny inferences can be used to:
- Improve genome annotation (by helping to identify mispredicted and/or missing genes)
- Help differentiate between orthologous and paralogous genes (by leveraging positional information)
- Perform evolutionary distance analyses (<i>e.g.</i> species sharing similar/identical genome reorganization events are likely more closely related than species that do not share these events).
</details>

<details open>
  <summary><b><i>Show/hide: Installation</i></b></summary>

## <b>Requirements</b>
- [DIAMOND](https://github.com/bbuchfink/diamond)
- [minimap2](https://github.com/lh3/minimap2)
- [MashMap3](https://github.com/marbl/MashMap)
- [Perl5](https://www.perl.org/)
- [PerlIO::gzip](https://metacpan.org/pod/PerlIO::gzip)
- [Python3](https://www.python.org/)
- [matplotlib](https://matplotlib.org/)
- [seaborn](https://seaborn.pydata.org/)
- [pandas](https://pandas.pydata.org/)
- [scipy](https://scipy.org/)
- [Circos](https://circos.ca/)

### <b>Downloading SYNY from GitHub</b>
```Bash
git clone https://github.com/PombertLab/SYNY.git
cd SYNY
export PATH=$PATH:$(pwd)
```

### <b>Installing SYNY</b>

SYNY dependencies can be installed automatically with Conda (does not require sudo privileges), with `setup_syny.pl` (requires sudo privileges), or manually.

<details open>
  <summary><b><i>Show/hide: Installing SYNY with Conda</i></b></summary>

#### <b>With Conda</b>
SYNY can be installed without sudo-elevated privileges by leveraging conda packages. The installation process was tested with Miniconda3 on Ubuntu-22.04.3 LTS and Fedora 40 Linux distributions running as virtual machines on Microsoft Windows Subsystem for Linux (WSL).

To install SYNY with conda:
```Bash
# Creating a conda environment
conda create -n syny  

# Activating the conda environment
conda activate syny

## Installing syny within the conda environment
conda install syny -c conda-forge -c bioconda 
```

To run SYNY within its conda environment:
```Bash
conda activate syny
(syny) username:~$ run_syny.pl -a *.gbff.gz -o output directory
```

If a conda package manager is needed, Miniconda can be installed without the need for sudo-elevated privileges. To install Miniconda (tested on Ubuntu/Fedora):
```Bash
## Setup variables
DOWN_DIR=$HOME/Downloads          ## Replace by desired download directory
TOOLS_DIR=$HOME/Tools             ## Replace by desired directory
CONDA_DIR=$TOOLS_DIR/miniconda3   ## Replace by desired subdirectory
CONFIG=~/.profile                 ## Desired configuration file

## Downloading miniconda
mkdir -p $DOWN_DIR $TOOLS_DIR
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -P $DOWN_DIR
chmod +x $DOWN_DIR/Miniconda3-latest-Linux-x86_64.sh

## Installing miniconda
$DOWN_DIR/Miniconda3-latest-Linux-x86_64.sh \
  -p $CONDA_DIR \
  -b

echo "export PATH=\$PATH:$CONDA_DIR/condabin" >> $CONFIG
source $CONFIG

## Deleting downloaded file after installation
rm $CONDA_DIR/Miniconda3-latest-Linux-x86_64.sh

## Initializing conda
conda init
bash --login

## Updating conda and installing the libmamba solver
conda update --yes -n base conda
conda install --yes -n base conda-libmamba-solver
conda config --set solver libmamba

## Disabling automatic activation of the base environment at login
conda config --set auto_activate_base false
```
</details>

<details open>
  <summary><b><i>Show/hide: Installing SYNY with Perl install script</i></b></summary>

#### <b>With Perl install script</b>
SYNY and its dependencies can be installed automatically with `setup_syny.pl`. This script will download and install [DIAMOND](https://github.com/bbuchfink/diamond), [minimap2](https://github.com/lh3/minimap2), [MashMap3](https://github.com/marbl/MashMap), [Circos](https://circos.ca/) together with the required dnf/apt/zypper packages (this script has been tested on Fedora, Ubuntu, Debian, Kali, and openSUSE Tumbleweed distributions). Note that using this script will require sudo privileges to install dnf/apt/zypper packages.

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
</details>



<details open>
  <summary><b><i>Show/hide: Installing SYNY manually</i></b></summary>

#### <b>Manual installation</b>
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

##### To install MashMap3:
```Bash
## On Ubuntu
sudo apt install build-essential
sudo apt install gsl-bin
sudo apt install libgsl-dev

## On Fedora
sudo dnf group install "C Development Tools and Libraries"
sudo dnf install cmake
sudo dnf install gsl gsl-devel

git clone https://github.com/marbl/MashMap.git
cd MashMap

cmake -H. -Bbuild -DCMAKE_BUILD_TYPE=Release
cmake --build build

mv ./build/bin ./bin

export PATH=$PATH:$(pwd)/bin ## MashMap3 install directory
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
  Text::Roman \
  Clone \
  Config::General \
  Font::TTF::Font \
  Getopt::ArgvFile \
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
</details>

</details>

<details open>
  <summary><b><i>Show/hide: About memory usage</i></b></summary>

### <b>About memory usage</b>
Performing pairwise alignments between genomes can quickly consume a large amount of RAM. The amount of memory required depends on the size of the genomes being aligned and on the amounts of repetitive elements present in their sequences. When working with genomes larger than 150 Mbp, we recommend using [MashMap3](https://github.com/marbl/MashMap) instead of [minimap2](https://github.com/lh3/minimap2) as the genome alignment tool. While MashMap3 does not product exact alignments <i>per se</i>, its results are congruent with the exact alignments produced by minimap2 (and as such constitute good approximations). For example, alignments between two repeat-laden 500 Mbp genomes ran within 5 GB of RAM with MashMap3 whereas the same alignments with minimap2 peaked around 100 GB of RAM despite the two tools producing similar outcomes. However, note that when using MashMap3, reducing its minimum percentage identity threshold below its 85% default settings will significantly increase its RAM usage. If running above the maximum amount of memory available (real + virtual memory), pairwise genome alignments will be terminated by Linux out-of-memory (OOM) killer, resulting in blank outputs for the corresponding alignments and plots.

</details>

<details open>
  <summary><b><i>Show/hide: Using SYNY</i></b></summary>

## <b>Using SYNY</b>
The SYNY pipeline can be run with [run_syny.pl](https://github.com/PombertLab/SYNY/blob/main/run_syny.pl), a master script that:
1. Extracts genome/protein sequences and annotation data from GenBank (.gbf/.gbff) flat files.
2. Performs round-robin pairwise genome alignments with [minimap2](https://github.com/lh3/minimap2) or [MashMap3](https://github.com/marbl/MashMap).
3. Performs round-robin [DIAMOND](https://github.com/bbuchfink/diamond) BLASTP homology searches, identifies conserved protein gene pairs, and reconstructs collinear clusters from these searches.
4. Generates dotplots, barplots, linemaps and [Circos](https://circos.ca/) plots highlighting collinear regions inferred from pairwise genome alignments and from shared protein cluster reconstructions.
5. Generates heatmaps summarizing the percentages of collinear bases found in pairwise genome alignments and the percentages of proteins found in collinear clusters between all genomes.

#### Why use two distinct approaches?
When working with well-annotated, closely related genomes, collinear regions inferred from pairwise genome alignments and shared protein cluster reconstructions should yield similar results.

However, when working with genomes featuring high levels of sequence divergence, pairwise genome alignments may struggle. In those instances, collinear regions inferred from protein cluster reconstructions should outperform those from genome alignments; hypervariable intergenic regions are not considered in protein cluster reconstructions and silent mutations/codon usage biases do not affect amino acid sequences. Amino acids also have a larger character state and lower back-mutation probabilities than nucleotide sequences, allowing for searches across wider evolutionary distances.

Conversely, when working with poorly annotated or unannotated genomes, collinear regions inferred from genome alignments should outperform those inferred from protein cluster reconstructions (if there is sufficient sequence similarity to perform pairwise alignments).

### Command line options
SYNY can be run from the master script as follows:<br>
```Bash
run_syny.pl \
  -a *.gbff \
  -o SYNY
```

The `run_syny.pl` command line options for can be entered directly from the command line and/or provided from one or more configuration file(s) containing one entry per line (see [commands.conf](https://github.com/PombertLab/SYNY/blob/main/Examples/commands.conf) for an example). For example:
```Bash
run_syny.pl \
  -a *.gbff \
  @commands.conf
```

Options for run_syny.pl are:
```
-h (--help)             Display all command line options
-t (--threads)          Number of threads to use [Default: 16]
-p (--pthreads)         Number of graphs to plot in parralel; defaults to --threads if unspecified
-a (--annot)            GenBank GBF/GBFF Annotation files (GZIP files are supported)
-o (--outdir)           Output directory [Default = SYNY]
-e (--evalue)           DIAMOND BLASTP evalue cutoff [Default = 1e-10]
-g (--gaps)             Allowable number of gaps between gene pairs [Default = 0]
--minsize               Minimum contig size (in bp) [Default: 1]
--include               Select contigs with names from input text file(s) (one name per line); i.e. exclude everything else
--ranges                Select contigs with subranges from input text file(s): name start end
--exclude               Exclude contigs with names matching the regular expression(s); e.g. --exclude '^AUX'
--aligner               Specify genome alignment tool: minimap or mashmap [Default: minimap]
--asm                   Specify minimap max divergence preset (--asm 5, 10 or 20) [Default: off]
--mpid                  Specify mashmap percentage identity [Default: 85]
--resume                Resume minimap/mashmap computations (skip completed alignments)
--min_asize             Filter out alignments/clusters smaller than integer value (e.g. --min_asize 5000)
--no_sec                Turn off minimap2 secondary alignments
--no_map                Skip minimap/mashmap pairwise genome alignments
--no_vcf                Skip minimap VCF file creation (files can be quite large)
--no_clus               Skip gene cluster reconstructions
--version               Display SYNY version

### Circos plots
-c (--circos)           Circos plot mode: pair (pairwise), cat (concatenated), all (cat + pair) [Default: pair]
--orientation           Karyotype orientation: normal, inverted or both [Default: normal]
--circos_prefix         Prefix for concatenated plots [Default: circos]
-r (--ref)              Reference to use for concatenated plots; uses first genome (alphabetically) if ommitted
-u (--unit)             Size unit (Kb or Mb) [Default: Mb]
--winsize               Sliding windows size (nucleotide biases) [Default: 10000]
--stepsize              Sliding windows step (nucleotide biases) [Default: 5000]
--labels                Contig label type: mixed (arabic + roman numbers), arabic, roman, or names [Default: mixed]
--label_size            Contig label size [Default: 36]
--label_font            Contig label font [Default: bold]
--custom_file           Load custom colors from file
--list_preset           List available custom color presets
--custom_preset         Use a custom color preset, e.g.: --custom_preset chloropicon
--max_ticks             Set max number of ticks [Default: 5000]
--max_ideograms         Set max number of ideograms [Default: 200]
--max_links             Set max number of links [Default: 75000]
--max_points_per_track  Set max number of points per track [Default: 75000]
--clusters              Color by cluster instead of contig/chromosome [Default: off]
--no_ntbiases           Turn off nucleotide biases and GC/AT skews subplots
--no_skews              Turn off GC / AT skews subplots
--no_cticks             Turn off ticks in Circos plots
--no_circos             Turn off Circos plots

### Barplots
-bh (--bheight)         Barplot figure height in inches [Default: 10.8]
-bw (--bwidth)          Barplot figure width in inches [Default: 19.2]
--bfsize                Barplot font size [Default: 8]
--palette               Barplot color palette [Default: Spectral]
--monobar               Use a monochrome barplot color instead: e.g. --monobar blue
--bclusters             Color clusters by alternating colors; colors are not related within/between contigs; [Default: off]
--bpmode                Barplot mode: pair (pairwise), cat (concatenated), all (cat + pair) [Default: pair]
--no_barplot            Turn off barplots

### Dotplots
-dh (--dheight)         Dotplot figure height in inches [Default: 10.8]
-dw (--dwidth)          Dotplot figure width in inches [Default: 19.2]
--dfsize                Dotplot font size [Default: 8]
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
--hfsize                Heatmap font size [Default: 8]
--hmpalette             Heatmap color palette [Default: winter_r]
--hmax                  Set maximum color bar value [Default: 100]
--hmin                  Set minimum color bar value [Default: 0]
--hauto                 Set color bar values automatically instead
--no_heatmap            Turn off heatmaps

### Linear maps
-lh (--lheight)         Linear map figure height in inches [Default: 5]
-lw (--lwidth)          Heatmap figure width in inches [Default: 20]
--lm_rpalette           Reference genome color palette [Default: tab20]
--lm_xpalette           Target genome color palette [Default: Blues]
--lmrotation            Contig name rotation [Default: 90]
--lfsize                Font size [Default: 8]
--no_linemap            Turn off linemaps
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
  - PAF: pairwise genome alignments in the corresponding format (minimap2/MashMap3)
  - MAF, ALN: pairwise genome alignments in the corresponding formats (minimap2)
  - VCF: Variant Call Format files from minimap2 alignments (can be turned off with --no_vcf)
  - METRICS (minimap2):
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
	  - Barplots (in PNG/SVG format) from minimap2/mashmap3 alignments (.mmap.)
	  - Barplots (in PNG/SVG format) from protein clusters found with SYNY (.gap_0., .gap_1., ...)
  - CIRCOS:
	  - Circos plots (in PNG/SVG format) from minimap2/mashmap3 alignments (.mmap.)
	  - Circos plots (in PNG/SVG format) from protein clusters found with SYNY (.gap_0., .gap_1., ...)
  - CIRCOS_DATA:
	  - Configuration files for Circos plots
  - DOTPLOTS:
  	- Dotplots (in PNG/SVG format) from minimap2/mashmap3 alignments (.mmap.)
  	- Dotplots (in PNG/SVG format) from protein clusters found with SYNY (.gap_0., .gap_1., ...)
  - HEATMAPS:
    - Heatmaps (in PNG/SVG format) summarizing the percentages of collinear bases in pairwise alignments (.mmap.)
    - Heatmaps (in PNG/SVG format) summarizing the percentages of proteins found in clusters (.gap_0., .gap_1., ...)
  - LINEMAPS:
    - Linear maps (in PNG/SVG format) from minimap2/mashmap3 alignments (.mmap.)
    - Linear maps (in PNG/SVG format) from protein clusters found with SYNY (.gap_0., .gap_1., ...)
- SEQUENCES:
  - GENOMES:
  	- FASTA files containing the sequences of the investigated genomes
  - PROTEINS:
    - Protein sequences for each species (.faa) in FASTA format
</details>

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

SYNY can be run from the command line with the `run_syny.pl` master script. In the command line below, no gap is allowed during gene cluster inferences, and Circos plots are produced in pairwise mode. The latter can be changed with the `--circos` command line switch (possible values are: `pair`, `concatenated`, `all`). Note that when comparing several genomes, concatenated plots can quickly become too dense for legibility. Producing concatenated plots can also significantly increase computation time.

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

To facilitate comparisons when working with large datasets, heatmaps summarizing the percentages of collinear bases in pairwise genome alignments (`.mmap.`) and the percentages of protein coding-genes found in collinear clusters between each pair of genomes (e.g. `.gap_0.`) are generated with matplotlib.

In the above example, small heatmaps with 25 datapoints (5 * 5 genomes) will be generated for the minimap2 pairwise alignments and for the protein clusters found for each gap value investigated. In these heatmaps, percentages between pairs of genomes will vary based on the total number of bases and proteins found in the query used: <i>i.e.</i> `(collinear bases / total bases in query) * 100` and `(proteins in clusters / total proteins in query) * 100`.

##### Example of a heatmap showing percentages of collinear protein-coding genes between investigated genomes:
<p align="left">
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/proteins_in_clusters.gap_0.heatmap.png">
</p>

Heatmap dimensions (default: 10 x 10) can be modified with the `--hheight` and `--hwidth` command line switches. Color bar min/max values (defaults: 0/100) can be changed with the `--hmin` and `--hmax` command line switches or calculated automatically with `--hauto`. The heatmap color palette (default: winter_r) can be changed with the `--hmpalette` command line switch (see this [URL](https://github.com/PombertLab/SYNY/blob/main/Images/python_color_palettes.png) for a list of color palettes). Color palettes available on the operating system can be listed and/or plotted with [check_mp_colors.py](https://github.com/PombertLab/SYNY/blob/main/Utils/check_mp_colors.py).

#### Circos plots

Unless the `--no_circos` command line switch is invoked, [Circos](https://circos.ca/) (pairwise and/or concatenated) will be generated from the protein clusters identified with SYNY (e.g. `.gap_0.`) and/or from the genome alignments computed with minimap2 or MashMap3 (`.mmap.`).

In the pairwise plots (`--circos pair`), genomes are plotted in pairs (query <i>vs.</i> subject) using the query as the reference. In the concatenated plots (`--circos cat`), all genomes are plotted together in a single figure, using the reference genome specified with the `--ref` command line switch. If omitted, the first genome encountered alphabetically will be used as the default reference. Both concatenated and pairwise plots can be generated with the `--circos all` command line switch.

Karyotypes can be plotted with Circos in normal and/or inverted orientation(s). In the inverted plots, the contigs/chromosomes from the genome(s) being compared to the reference are plotted in reverse, from last to first. By default, the Circos plots are plotted in normal orientation. This behavior can be modified with the `--orientation` command line switch (possible values are `normal`, `inverted` and `both`).

##### Example of an image generated with Circos and SYNY from shared proteins clusters:
<p align="left">
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/WM276_vs_JEC21.gap_0.normal.png">
</p>

In this figure, nucleotides biases and GC/AT skews are plotted in the concentric rings. These subplots can be skipped with the `--no_ntbiases` command line switch. GC/AT skews can be further turned off independently from nucleotide biases subplots with the `--no_skews` options.

From from outer to inner, the concentric rings are:

- AT and GC nucleotide biases (grey and red lines)
- GT and AC nucleotide biases (blue and green lines)
- GA and CT nucleotide biases (purple and yellow lines)
- GC skews (histogram: blue [positive], red [negative] values)
- AT skews (histogram: orange [positive], purple [negative] values)


Syntenic blocks identified by SYNY are indicated by ribbons. These ribbons are color-coded based on the chromosomes/contigs present in the reference genome used. If desired, color coding can be set by cluster instead with `--clusters`. This option is useful when working with prokaryotes featuring a single chromosome, so that the ribbons are not all of the same color.

By default, contigs will be labelled by mixed numerals (roman + arabic) in the Circos plots. Contigs from the reference (query) will be labelled by `roman` numerals, others by `arabic` numerals. This behavior can be changed with the `--labels` command line option. Possible values are `mixed`, `roman`, `arabic` and `names`. Using `--labels roman` or `--labels arabic` will set all label numbers to the corresponding format, whereas using `--labels names` will label contigs by their names instead. Label sizes (default: 36) and fonts can be further adjusted with the `--label_size` and `--label_font` command line options. Possible fonts are: `light`, 
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
  --no_skews \
  --labels names \
  --label_size 20 \
  --label_font semibold
```

##### Example of the same Circos plot, this time labelled by contig names intead:
<p align="left">
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/WM276_vs_JEC21.gap_0.normal.names.png">
</p>

All data and configuration files for the Circos plots are in the `PLOTS/CIRCOS_DATA/` subdirectory. If desired, chromosome labels can be further adjusted manually by editing the `LABEL` column in the corresponding Circos [karyotype](https://circos.ca/documentation/tutorials/ideograms/karyotypes/) file(s). Karyotype files created by SYNY (in `PLOTS/CIRCOS_DATA/` ) ends with `*.karyotype`.

```Bash
## To see the first 5 lines from a karyotype file:
head -n 5 concatenated.normal.karyotype

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
Chromosome maps (aka barplots) highlighting collinear genome segments will be generated with [paf_to_barplot.py](https://github.com/PombertLab/SYNY/blob/main/Plots/paf_to_barplot.py) and [matplotlib](https://matplotlib.org/) from minimap2/mashmap3 pairwise genome alignments (`.mmap.`) and from protein clusters identified with SYNY (e.g. `.gap_0.`).

##### Example of a barplot generated from minimap2 PAF files using defaults settings:
<p align="left">
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/WM276_vs_JEC21.mmap.barplot.19.2x10.8.Spectral.png">
</p>

In these plots, collinear regions found between the compared genomes are highlighted by colored rectangles. By default, these rectangles are color-coded based on the contigs/chromosomes of the query. The above barplot image was generated using the Spectral color palette from [seaborn](https://seaborn.pydata.org/tutorial/color_palettes.html) (set as default in SYNY). This palette can be replaced using the `--palette` command line switch; <i>e.g.</i> `--palette husl` (see this [URL](https://github.com/PombertLab/SYNY/blob/main/Images/python_color_palettes.png) for a list of color palettes). Color palettes available on the operating system can be listed and/or plotted with [check_mp_colors.py](https://github.com/PombertLab/SYNY/blob/main/Utils/check_mp_colors.py).

If desired, the barplots can instead be generated using a single monochromatic color with the `--monobar` command line switch; <i>e.g.</i> `--monobar red`. If this option is selected, the color-coding legend will be omitted from the plot.

By default, the barplots are formatted for a widescreen (landscape) output (width/height ratio: 19.2/10.8). This ratio can be adjusted with the `--bheight` and `--bwidth` command line switches.

By default, SYNY generates pairwise barplots. If desired, concatenated barplots can also be generated. In these concatenated barplots, the contigs from all species queried against a species (the reference) are plotted together. This behaviour can be changed with the `--bpmode` option. Possible values are `pair` (pairwise; set as default), `cat` (concatenated), and `all` (cat + pair).

#### Linemaps
Linear chromosome maps (aka linemaps) highlighting collinear genome segments will be generated with [linear_maps.py](https://github.com/PombertLab/SYNY/blob/main/Plots/linear_maps.py) and [matplotlib](https://matplotlib.org/) from pairwise genome alignments (`.mmap.`) and from protein clusters identified with SYNY (e.g. `.gap_0.`). In these plots, collinear regions found between the compared genomes are highlighted by colored polygons. By default, these polygons are color-coded based on the contigs/chromosomes of the reference. The above linemap image was generated using the tab20/Blues color palettes for the reference/target genomes (both set as defaults in SYNY). These palettes can be replaced using the `--lm_rpalette` and `--lm_xpalette` command line switches; <i>e.g.</i> `--lm_xpalette Grays` (see [URL](https://github.com/PombertLab/SYNY/blob/main/Images/python_color_palettes.png) for a list of color palettes). By default, the linemaps are formatted using a 20:5 width:height ratio. This ratio can be adjusted with the `--lheight` and `--lwidth` command line switches.

##### Example of a linemap generated from gene clusters (no gaps allowed) using defaults settings:
<p align="left">
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/WM276_vs_JEC21.gap_0.linemap.20x5.tab20.png">
</p>

#### Dotplots

Unless the `--no_dotplot` command line switch is invoked, dotplots will be generated with [paf_to_dotplot.py](https://github.com/PombertLab/SYNY/blob/main/Plots/paf_to_dotplot.py) and [matplotlib](https://matplotlib.org/) from pairwise genome alignments (`.mmap.`) and from protein clusters identified with SYNY, (e.g. `.gap_0.`).

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

The default monochromatic color (blue) can be changed with the `--color` option (<i>e.g.</i> `--color red`). If desired, dotplots can instead be color-coded based on the query contigs/chromosomes with the `--dotpalette` option; <i>e.g.</i> `--dotpalette husl` (see this [URL](https://github.com/PombertLab/SYNY/blob/main/Images/python_color_palettes.png) for a list of color palettes). Color palettes available on the operating system can be listed and/or plotted with [check_mp_colors.py](https://github.com/PombertLab/SYNY/blob/main/Utils/check_mp_colors.py).

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

PAF metrics cannot be calculated for [MashMap3](https://github.com/marbl/MashMap) aligments since the data required is missing from its PAF output files. MashMap3 does not calculate alignments explicitly (as per its intructional manual) but runs in a much smaller memory footprint than minimap2 when using its default percentage identity (`--mpid 85`). As such, it constitutes an interesting alternative when running into memory constraints. However, note that lowering the MashMap3 default threshold will significantly increase its memory usage.

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

Here, while the pairwise alignments and gene cluster collinear inferences produce similar results within species (not shown), they differ substantially between species. For example, below are two barplots between <i>Encephalitozoon hellem</i> and <i>Ordospora colligata</i>. While the results are congruent, gene clusters inferences `.gap_0.` perform better than those based on pairwise alignments (`.mmap.`), the latter of which struggle due to the high levels of sequence divergence involved. These differences will also be reflected in the dotplots, heatmaps and Circos plots generated (not shown)

When working with genomes exhibiting a high level of divergence, both minimap2 and MashMap3 may fail to produce decent alignments. If that happens, we recommend using MashMap3 with a minimum percentage identity set lower than its default 85% threshold (e.g. `--mpid 70`). Note that doing so will increase its RAM usage and may introduce artefacts/noise in the plots produced.

Barplots (gene clusters): 
<p align="left">
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/Eh50604_vs_OcOC4.gap_0.barplot.19.2x10.8.Spectral.png">
</p>

Barplots (genome alignments; minimap2): 
<p align="left">
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/Eh50604_vs_OcOC4.mmap.barplot.19.2x10.8.Spectral.png">
</p>

Barplots (genome alignments; mashmap; --mpid 85): 
<p align="left">
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/Eh50604_vs_OcOC4.mmap85.barplot.19.2x10.8.Spectral.png">
</p>

Barplots (genome alignments; mashmap; --mpid 70): 
<p align="left">
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/Eh50604_vs_OcOC4.mmap70.barplot.19.2x10.8.Spectral.png">
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
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/encephalitozoon.gap_0.normal.png">
</p>

In the above [Circos](https://circos.ca/) plot, links between the reference and the queried genomes are color-coded based on the reference karyotype; links between the queried genomes are shown in light gray. In these plots, each chromosome/contig is plotted as a distinct [ideogram](https://circos.ca/documentation/tutorials/ideograms/ideograms/).

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
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/encephalitozoon_cc.gap_0.normal.png">
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
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/encephalitozoon_blues.gap_0.normal.png">
</p>

### <b>File conversion</b>
##### <i>Fasta to GBFF</i>
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

##### <i>Fasta + GFF3 to GBFF</i>
A simple Fasta + GFF3 to GBFF converter ([gff3_to_gbff.pl](https://github.com/PombertLab/SYNY/blob/main/Utils/gff3_to_gbff.pl)) is available in the `Utils/` subdirectory. This tool was tested on NCBI GFF3 files and expects the GFF3 file(s) to include gene/mRNA/exon/CDS entries in the `type` column and the `ID` and `Parent` tags in the attributes column. It also expects the corresponding Fasta and GFF3 files to share the same prefixes (<i>e.g.</i> genome_1.fasta / genome_1.gff). The GBFF files created by [gff3_to_gbff.pl](https://github.com/PombertLab/SYNY/blob/main/Utils/gff3_to_gbff.pl) were designed to work with SYNY but do not adhere exactly to the GBFF format and may not work for other purposes.

To convert FASTA + GFF3 file(s) to pseudo-GBFF files (compressed in gzip format) using the standard genetic code:
```Bash
gff3_to_gbff.pl \
  --fasta FASTA/* \
  --gff3 GFF/* \
  --outdir GBFF \
  --gzip \
  --gcode 1 \
  --verbose
```

Options for `gff3_to_gbff.pl` are:
```
-f (--fasta)    FASTA file(s) to convert (gziped files are supported)
-g (--gff3)     GFF3 files to convert (gziped files are supported)
-o (--outdir)   Output directory [Default: GBFF]
-z (--gzip)     Compress the GBFF output files
-c (--gcode)    NCBI genetic code [Default: 1]
                1  - The Standard Code
                2  - The Vertebrate Mitochondrial Code
                3  - The Yeast Mitochondrial Code
                4  - The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
                11 - The Bacterial, Archaeal and Plant Plastid Code
                NOTE - For complete list; see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
-v (--verbose)  Add verbosity
--version       Show script version
```
When compressing the GBFF output files, gff3_to_gbff.pl will use pigz if available, otherwise it will default to gzip

##### <i>JGI GFF to NCBI GFF3 format</i>
A simple JGI GFF to NCBI GFF3 converter ([jgi_to_ncbi_gff.pl](https://github.com/PombertLab/SYNY/blob/main/Utils/jgi_to_ncbi_gff.pl)) is available in the `Utils/` subdirectory. This tool converts JGI GFF files to pseudo-NCBI GFF3 files compatible with [gff3_to_gbff.pl](https://github.com/PombertLab/SYNY/blob/main/Utils/gff3_to_gbff.pl).

To convert JGI GFF to pseudo-GFF3 files:
```Bash
jgi_to_ncbi_gff \
  --gff3 jgi.gff \
  --out ncbi.gff3
```

Options for `jgi_to_ncbi_gff.pl` are:
```
-g (--gff3)     JGI GFF files to convert
-o (--out)      Output file
--version       Show script version
```

##### <i>Reordering/reorienting contigs</i>
When working with newly assembled genomes, contigs are sometimes found out-of-order and/or on opposite strands when compared to reference genomes. As such, reordering/reorienting contigs based on their reference genome counterparts is often useful to facilitate comparisons. A simple python script [orient_fastas_to_reference.py](https://github.com/PombertLab/SYNY/blob/main/Utils/orient_fastas_to_reference.py) is available in the `Utils/` subdirectory to help with this task. In a nutshell, `orient_fastas_to_reference.py` performs BLASTN homology searches between genomes (queries) and a reference genome (subject) by leveraging [NCBI BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/), then reorder/reorient contigs from the queries according to the results of these homology searches.

To reorder/reorient contigs with `orient_fastas_to_reference.py`:
```Bash
orient_fastas_to_reference.py \
  --fasta *.fasta \
  --ref reference.fasta \
  --outdir FASTA_oriented \
  --min_pident 75 \
  --min_palign 20 \
  --max_overlp
```

Options for `orient_fastas_to_reference.py` are:
```
-f (--fasta)         FASTA files to reorder/reorient
-r (--ref)           Reference genome assembly
-o (--outdir)        Output directory [Default:'oriented_fastas']
-i (--min_pident)    Minimum percent identity to assign segment to reference [Default: 85%]
-a (--min_palign)    Minimum percent of the contig participating in alignment to assign segment to reference [Default: 20%]
-x (--max_overlp)    Maximum percent of alignment allowed to overlap a previous alignment to assign segment to reference [Default: 5%]
--version            Show script version
```

A subdirectory will be created foreach FASTA file queried. Inside each subdirectory, BLASTN results will be located in `results.blastn.6`. Contigs with or without matches against the reference genomes will be located in the `.oriented.fasta` and `unmatched.fasta` files, respectively, and summaries will be found in `all.map`. Circos karyotype and links files will also be generated automatically.

### Example 3 - Subsets
Because large and/or fragmented genomes can be difficult to visualize due to plot density, SYNY also includes options to look at subsets of genomes. Contigs to be investigated/plotted can be specified from a single text file with the `--include` command line switch, while portions of contigs can also be specified in a tab-delimited text file with the `--ranges` command line switch. Alternatively, contigs can also be excluded using a regular expression with the `--exclude`command line switch. The latter option is useful when dealing with accession numbers containing a mixture of complete chromosomes and partial contigs, the latter often indicated with distinctive names. These partial contigs are often small, numerous, and tend to clutter plots. As such, removing them is often desirable.

Below are example to plots portions from two Arabidopsis genomes (<i>A. thaliana</i> and <i>A. arenosa</i>). These genomes can be downloaded automatically using the Arabidopsis.sh shell script (from the `Examples/` directory).

```Bash
## Downloading example data from NCBI using the Arabidopsis.sh shell script 

DATA=~/DATA               ## Replace ~/DATA by desired annotation data directory
Arabidopsis.sh $DATA
```

To plot the full dataset:

```Bash
##### Runtime (Intel i5-12500H mobile CPU)
# Total: 6 min; Circos: 4 min)

DATA=~/DATA                   ## Replace by annotation data directory
SYNY=~/SYNY_RESULTS/ARAB_ALL  ## Replace by desired output directory

run_syny.pl \
  --threads 16 \
  --annot $DATA/*.gbff.gz \
  --aligner mashmap \
  --outdir $SYNY \
  --g 0 1 5
```

<p align="left">
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/TAIR10_vs_AARE701.gap_0.barplot.19.2x10.8.Spectral.png">
</p>

In the corresponding plots, two small contigs (NC_000932 and NC_037304) are observed. These are organelle genomes (chloroplast + mitochondrion) sometimes present together with nuclear genomes in accesssion numbers.

To plot only the chromosomes originating from the nuclear genomes, we can specify them in a simple text file, then run SYNY with the `--include` command line switch:

```Bash
# Content of contigs.txt
NC_003070
NC_003071
NC_003074
NC_003075
NC_003076
LR999451
LR999452
LR999453
LR999454
LR999455
LR999456
LR999457
LR999458
```

```Bash
##### Runtime (Intel i5-12500H mobile CPU)
# Total: 6 min; Circos: 4 min)

DATA=~/DATA                   ## Replace by annotation data directory
SYNY=~/SYNY_RESULTS/ARAB_INC  ## Replace by desired output directory

run_syny.pl \
  --threads 16 \
  --annot $DATA/*.gbff.gz \
  --aligner mashmap \
  --outdir $SYNY \
  --g 0 1 5 \
  --include contigs.txt
```

If desired we can also plot only segments from specified contigs/chromosomes. These segments can be specified in a simple (tab/space)-delimited text file, then run SYNY with the `--ranges` command line switch:

```Bash
# Content of subranges.txt
# Contig    start      end
NC_003070   1          1000000
NC_003070   1500000    2500000
NC_003070   3000000    4000000
NC_003071   1          19698289
LR999451    1          25224288
LR999452    1          14316399
LR999453    1          21694789
```

```Bash
##### Runtime (Intel i5-12500H mobile CPU)
# Total: 1 min; Circos: 30 sec)

DATA=~/DATA                   ## Replace by annotation data directory
SYNY=~/SYNY_RESULTS/ARAB_SUB  ## Replace by desired output directory

run_syny.pl \
  --threads 16 \
  --annot $DATA/*.gbff.gz \
  --aligner mashmap \
  --outdir $SYNY \
  --g 0 1 5 \
  --ranges subranges.txt
```

In the above example, numbers will be appended to contig names whenever approriate to differentiate between multiple segments, i.e. NC_003070_1, NC_003070_2 and NC_003070_3 (see below).

<p align="left">
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/TAIR10_vs_AARE701.mmap.barplot.19.2x10.8.Spectral.png">
</p>

<details open>
  <summary><b><i>Show/hide: Funding and acknowledgments</i></b></summary>

## Funding and acknowledgments
This work was supported in part by the National Institute of Allergy and Infectious Diseases of the National Institutes of Health (award number R15AI128627) to Jean-Francois Pombert. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.
</details>

<details open>
  <summary><b><i>Show/hide: How to cite</i></b></summary>

## How to cite
##### A short paper describing the SYNY pipeline is available as a preprint on bioRxiv:

[SYNY: a pipeline to investigate and visualize collinearity between genomes](https://www.biorxiv.org/content/10.1101/2024.05.09.593317v1). Julian AT, Pombert JF. <b>bioRxiv</b>, 2024.05.09.593317. DOI: 10.1101/2024.05.09.593317.

##### If you use SYNY, please also cite the tool(s) used for genome alignments and/or protein sequence homology, as needed:

[Minimap2: pairwise alignment for nucleotide sequences](https://pubmed.ncbi.nlm.nih.gov/29750242/). Li H. <b>Bioinformatics.</b> 2018 Sep 15;34(18):3094-3100. doi: 10.1093/bioinformatics/bty191.

[Minmers are a generalization of minimizers that enable unbiased local Jaccard estimation](https://www.biorxiv.org/content/10.1101/2023.05.16.540882v1). Kille B, Garrison E, Treangen TJ, Phillippy AM. <b>Bioinformatics</b>. 2023 Sep 2;39(9):btad512. PMID: 37603771 PMCID: PMC10505501. doi: 10.1093/bioinformatics/btad512.

[Sensitive protein alignments at tree-of-life scale using DIAMOND](https://www.nature.com/articles/s41592-021-01101-x). Buchfink B, Reuter K, Drost HG. <b>Nature Methods.</b> 18, 366–368 (2021). doi: 10.1038/s41592-021-01101-x
</details>

<details open>
  <summary><b><i>Show/hide: References</i></b></summary>

## <b>References</b>
[Sensitive protein alignments at tree-of-life scale using DIAMOND](https://www.nature.com/articles/s41592-021-01101-x). Buchfink B, Reuter K, Drost HG. <b>Nature Methods.</b> 18, 366–368 (2021). doi: 10.1038/s41592-021-01101-x

[Minimap2: pairwise alignment for nucleotide sequences](https://pubmed.ncbi.nlm.nih.gov/29750242/). Li H. <b>Bioinformatics.</b> 2018 Sep 15;34(18):3094-3100. doi: 10.1093/bioinformatics/bty191.

[Minmers are a generalization of minimizers that enable unbiased local Jaccard estimation](https://www.biorxiv.org/content/10.1101/2023.05.16.540882v1). Kille B, Garrison E, Treangen TJ, Phillippy AM. <b>Bioinformatics</b>. 2023 Sep 2;39(9):btad512. PMID: 37603771 PMCID: PMC10505501. doi: 10.1093/bioinformatics/btad512.

[Circos: an information aesthetic for comparative genomics](https://pubmed.ncbi.nlm.nih.gov/19541911/). Krzywinski M, Schein J, Birol I, Connors J, Gascoyne R, Horsman D, Jones SJ, Marra MA. <b>Genome Res.</b> 2009 Sep;19(9):1639-45. doi: 10.1101/gr.092759.109.

[Matplotlib: A 2D graphics environment](https://doi.org/10.1109/MCSE.2007.55). Hunter, JD. <b>Computing in Science & Engineering.</b> 2007 9(3):90-95. doi: 10.1109/MCSE.2007.55.

[Seaborn: statistical data visualization](https://doi.org/10.21105/joss.03021). Waskom, ML. <b>Journal of Open Source Software.</b> 2021 6(60): 3021. doi: 10.21105/joss.03021

[Data structures for statistical computing in Python](https://doi.org/10.25080/Majora-92bf1922-00a). McKinney W. <b>Proceedings of the 9th Python in Science Conference.</b> 2010:56-61. doi: 10.25080/Majora-92bf1922-00a.

[SciPy 1.0: fundamental algorithms for scientific computing in Python.](https://pubmed.ncbi.nlm.nih.gov/32015543/). Virtanen P, Gommers R, Oliphant TE, Haberland M, Reddy T, Cournapeau D, Burovski E, Peterson P, Weckesser W, Bright J, van der Walt SJ, Brett M, Wilson J, Millman KJ, Mayorov N, Nelson ARJ, Jones E, Kern R, Larson E, Carey CJ, Polat İ, Feng Y, Moore EW, VanderPlas J, Laxalde D, Perktold J, Cimrman R, Henriksen I, Quintero EA, Harris CR, Archibald AM, Ribeiro AH, Pedregosa F, van Mulbregt P; SciPy 1.0 Contributors. <b>Nat Methods.</b> 2020 Mar;17(3):261-272. doi: 10.1038/s41592-019-0686-2. Epub 2020 Feb 3. PMID: 32015543; PMCID: PMC7056644.

[The genome of the basidiomycetous yeast and human pathogen <i>Cryptococcus neoformans</i>](https://pubmed.ncbi.nlm.nih.gov/15653466/). Loftus BJ, Fung E, Roncaglia P, Rowley D, Amedeo P, Bruno D, Vamathevan J, Miranda M, Anderson IJ, Fraser JA, Allen JE, Bosdet IE, Brent MR, Chiu R, Doering TL, Donlin MJ, D'Souza CA, Fox DS, Grinberg V, Fu J, Fukushima M, Haas BJ, Huang JC, Janbon G, Jones SJ, Koo HL, Krzywinski MI, Kwon-Chung JK, Lengeler KB, Maiti R, Marra MA, Marra RE, Mathewson CA, Mitchell TG, Pertea M, Riggs FR, Salzberg SL, Schein JE, Shvartsbeyn A, Shin H, Shumway M, Specht CA, Suh BB, Tenney A, Utterback TR, Wickes BL, Wortman JR, Wye NH, Kronstad JW, Lodge JK, Heitman J, Davis RW, Fraser CM, Hyman RW. <b>Science.</b> 2005 Feb 25;307(5713):1321-4. doi: 10.1126/science.1103773. Epub 2005 Jan 13. PMID: 15653466; PMCID: PMC3520129.

[Genome variation in <i>Cryptococcus gattii</i>, an emerging pathogen of immunocompetent hosts](https://pubmed.ncbi.nlm.nih.gov/21304167/). D'Souza CA, Kronstad JW, Taylor G, Warren R, Yuen M, Hu G, Jung WH, Sham A, Kidd SE, Tangen K, Lee N, Zeilmaker T, Sawkins J, McVicker G, Shah S, Gnerre S, Griggs A, Zeng Q, Bartlett K, Li W, Wang X, Heitman J, Stajich JE, Fraser JA, Meyer W, Carter D, Schein J, Krzywinski M, Kwon-Chung KJ, Varma A, Wang J, Brunham R, Fyfe M, Ouellette BF, Siddiqui A, Marra M, Jones S, Holt R, Birren BW, Galagan JE, Cuomo CA. <b>mBio</b>. 2011 Feb 8;2(1):e00342-10. doi: 10.1128/mBio.00342-10. PMID: 21304167; PMCID: PMC3037005.

[A New Lineage of <i>Cryptococcus gattii</i> (VGV) Discovered in the Central Zambezian Miombo Woodlands](https://pubmed.ncbi.nlm.nih.gov/31719178/). Farrer RA, Chang M, Davis MJ, van Dorp L, Yang DH, Shea T, Sewell TR, Meyer W, Balloux F, Edwards HM, Chanda D, Kwenda G, Vanhove M, Chang YC, Cuomo CA, Fisher MC, Kwon-Chung KJ. <b>mBio</b>. 2019 Nov 12;10(6):e02306-19. doi: 10.1128/mBio.02306-19. PMID: 31719178; PMCID: PMC6851281.

[Telomere-to-Telomere genome assemblies of human-infecting <i>Encephalitozoon</i> species](https://pubmed.ncbi.nlm.nih.gov/37142951/). Mascarenhas Dos Santos AC, Julian AT, Liang P, Juárez O, Pombert JF. <b>BMC Genomics</b>. 2023 May 4;24(1):237. doi: 10.1186/s12864-023-09331-3. PMID: 37142951; PMCID: PMC10158259.

[A new microsporidian parasite, <i>Ordospora pajunii</i> sp. nov (Ordosporidae), of <i>Daphnia longispina</i> highlights the value of genomic data for delineating species boundaries](https://pubmed.ncbi.nlm.nih.gov/35279911/). de Albuquerque NRM, Haag KL, Fields PD, Cabalzar A, Ben-Ami F, Pombert JF, Ebert D. <b>J Eukaryot Microbiol</b>. 2022 May;69(3):e12902. doi: 10.1111/jeu.12902. Epub 2022 Mar 28. PMID: 35279911.

[The <i>Ordospora colligata</i> genome: Evolution of extreme reduction in microsporidia and host-to-parasite horizontal gene transfer](https://pubmed.ncbi.nlm.nih.gov/25587016/). Pombert JF, Haag KL, Beidas S, Ebert D, Keeling PJ. <b>mBio</b>. 2015 Jan 13;6(1):e02400-14. doi: 10.1128/mBio.02400-14. PMID: 25587016; PMCID: PMC4313915.
</details>
