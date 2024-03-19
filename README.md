## <b>Synopsis</b>

The SYNY pipeline investigates gene colinearity (synteny) between genomes by reconstructing clusters from conserved pairs of protein-coding genes identified from [DIAMOND](https://github.com/bbuchfink/diamond) homology searches. It also infers colinearity from pairwise genome alignments with [minimap2](https://github.com/lh3/minimap2).

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
    * [Example 2: <i>Encephalitozoon</i>](#Example-2---Encephalitozoon)
    * [Example 3: <i>Encephalitozoon</i> (with custom colors)](#Example-3---Encephalitozoon-with-custom-colors)
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
- [minimap2](https://github.com/lh3/minimap2)
- [Perl5](https://www.perl.org/)
- [PerlIO::gzip](https://metacpan.org/pod/PerlIO::gzip)
- [Python3](https://www.python.org/)
- [matplotlib](https://matplotlib.org/)
- [seaborn](https://seaborn.pydata.org/)

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
## On Ubuntu:
sudo apt install libperlio-gzip-perl

## On Fedora:
sudo dnf install perl-PerlIO-gzip
```


##### To install matplolib/seaborn:
```Bash
## On Ubuntu:
sudo apt install python3-matplotlib
sudo apt install python3-seaborn

## On Fedora:
sudo dnf install python3-matplotlib
sudo dnf install python3-seaborn

## Or via pip (Ubuntu/Fedora):
pip install matplotlib
pip install seaborn
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

printf "\nexport PATH=$PATH:$(pwd)" >> ~/.profile      ## Ubuntu
printf "\nexport PATH=$PATH:$(pwd)" >> ~/.bash_profile ## Fedora
```

## <b>Using SYNY</b>

The SYNY pipeline can be run with [run_syny.pl](https://github.com/PombertLab/SYNY/blob/main/run_syny.pl), a master script that:
1. Extracts genome and protein sequences from GenBank (.gbf/.gbff) annotation files.
2. Performs round-robin pairwize genome alignments with [minimap2](https://github.com/lh3/minimap2) and plots them with [matplotlib](https://matplotlib.org/).
3. Performs round-robin [DIAMOND](https://github.com/bbuchfink/diamond) BLASTP homology searches, identifies conserved protein gene pairs, and reconstructs colinear clusters from these searches.
4. Generates [Circos](https://circos.ca/) plots highlighting colinear regions inferred from pairwise genome alignments and from shared protein cluster reconstructions.

#### Why use two distinct approaches?

When working with well-annotated, closely related genomes, colinear regions inferred from pairwise genome alignments and shared protein cluster reconstructions should yield similar results.

However, when working with genomes featuring high levels of sequence divergence, pairwise genome alignments may struggle. In those instances, colinear regions inferred from protein cluster reconstructions should outperform those from genome alignments; silent mutations/codon usage biases do not affect amino acid sequences and hypervariable intergenic regions are not considered in protein cluster reconstructions.

Conversely, when working with poorly annotated or unannotated genomes, colinear regions inferred from genome alignments should outperform those inferred from protein cluster reconstructions (if there is sufficient sequence similarity to perform pairwise alignments).

### Command line options

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
-a (--annot)    GenBank GBF/GBFF Annotation files (GZIP files are supported)
-e (--evalue)   BLAST evalue cutoff [Default = 1e-10]
-g (--gaps)     Allowable number of gaps between pairs [Default = 0]
-o (--outdir)   Output directory [Default = SYNY]
--asm           Specify minimap2 max divergence preset (--asm 5, 10 or 20) [Default: off]
--no_map        Skip minimap2 pairwise genome alignments

OPTIONS (PLOTS):
### Circos
-c (--circos)   Generate Circos plots automatically - http://circos.ca/
-circos_prefix  Desired Circos plot prefix [Default: circos]
-r (--ref)      Genome to use as reference (defaults to first one alphabetically if none provided)
-u (--unit)     Size unit (Kb or Mb) [Default: Mb]
--winsize       Sliding windows size (nucleotide biases) [Default: 10000]
--stepsize      Sliding windows step (nucleotide biases) [Default: 5000]
-custom_file    Load custom colors from file
-list_preset    List available custom color presets
-custom_preset  Use a custom color preset; e.g.
                # chloropicon - 20 colors - Lemieux et al. (2019) https://pubmed.ncbi.nlm.nih.gov/31492891/
                # encephalitozoon - 11 colors - Pombert et al. (2012) https://pubmed.ncbi.nlm.nih.gov/22802648/
-max_ticks              Set max number of ticks [Default: 5000]
-max_ideograms          Set max number of ideograms [Default: 200]
-max_links              Set max number of links [Default: 25000]
-max_points_per_track   Set max number of points per track [Default: 75000]

### Barplots/Dotplots
-h (--height)   Figure height in inches [Default: 10.8]
-w (--width)    Figure width in inches [Default: 19.2]
-m (--multi)    Axes units multiplier (for dotplots) [Default: 1e5]
--palette       Color palette (for barplots) [Default: Spectral]
--monobar       Use a monochrome color for barplots instead of color palette: e.g. --monobar blue
--color         Scatter plot color (for dotplots) [Default: blue]
--dotpalette    Use a color palette instead of a monochrome color for dotplots: e.g. --dotpalette inferno
--noticks       Turn off ticks on x and y axes
--wdis          Horizontal distance (width) between subplots [Default: 0.05]
--hdis          Vertical distance (height) between subplots [Default: 0.1]
--no_dotplot    Skip dotplot creation
```
The output directory will be structured as follows: 
```Bash
ls -lah SYNY/

drwxr-xr-x 11 jpombert jpombert 4.0K Mar  8 14:46 .
drwx------ 23 jpombert jpombert 4.0K Mar  8 14:45 ..
drwxr-xr-x  5 jpombert jpombert 4.0K Mar  8 14:45 ALIGNMENTS
drwxr-xr-x  3 jpombert jpombert 4.0K Mar  8 14:45 BARPLOTS
drwxr-xr-x  9 jpombert jpombert 4.0K Mar  8 14:45 CIRCOS
drwxr-xr-x  2 jpombert jpombert 4.0K Mar  8 14:45 CONSERVED
drwxr-xr-x  3 jpombert jpombert 4.0K Mar  8 14:45 DIAMOND
drwxr-xr-x  3 jpombert jpombert 4.0K Mar  8 14:45 DOTPLOTS
drwxr-xr-x  2 jpombert jpombert 4.0K Mar  8 14:45 GENOME
drwxr-xr-x  2 jpombert jpombert 4.0K Mar  8 14:45 LISTS
drwxr-xr-x  2 jpombert jpombert 4.0K Mar  8 14:45 PROT_SEQ
drwxr-xr-x  2 jpombert jpombert 4.0K Mar  8 14:45 SHARED
drwxr-xr-x  5 jpombert jpombert 4.0K Mar  8 14:45 SYNTENY
-rw-r--r--  1 jpombert jpombert  182 Mar  8 14:45 error.log
-rw-r--r--  1 jpombert jpombert  397 Mar  8 14:46 syny.log
-rw-r--r--  1 jpombert jpombert 2.5M Mar  8 14:46 circos.paf.inverted.png
-rw-r--r--  1 jpombert jpombert 3.4M Mar  8 14:46 circos.paf.inverted.svg
-rw-r--r--  1 jpombert jpombert 2.4M Mar  8 14:46 circos.paf.normal.png
-rw-r--r--  1 jpombert jpombert 3.4M Mar  8 14:46 circos.paf.normal.svg
-rw-r--r--  1 jpombert jpombert 1.2M Mar  8 14:46 circos.syny.inverted.png
-rw-r--r--  1 jpombert jpombert 3.3M Mar  8 14:46 circos.syny.inverted.svg
-rw-r--r--  1 jpombert jpombert 1.1M Mar  8 14:45 circos.syny.normal.png
-rw-r--r--  1 jpombert jpombert 3.3M Mar  8 14:45 circos.syny.normal.svg
```

In the above, Circos colinearity plots inferred from pairwise genome alignments with [minimap2](https://github.com/lh3/minimap2) are indicated with the `.paf.` tag. Colinearity plots inferred from shared protein clusters identified with [SYNY](https://github.com/PombertLab/SYNY) are indicated with the `.syny.` tag.

The contents of the subdirectories are:
- ALIGNMENTS:
	- Pairwise minimap2 genome alignments in MAF, PAF and ALN formats
- BARPLOTS:
	- Barplots (in PNG/SVG format) generated from the minimap2 PAF alignments
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
- DOTPLOTS:
	- Dotplots (in PNG format) generated from the minimap2 PAF alignments
- GENOME:
	- FASTA files containing the sequences of the investigated genomes
- LISTS:
	- Lists of protein coding genes with location details (.list)
- PROT_SEQ:
	- Protein sequences for each species (.faa) in FASTA format
- SHARED:
	- Lists of all proteins and their top homologs (if any) in other species (.shared.tsv)
	- Lists of proteins that are unique to each species (.uniques.tsv)
- SYNTENY:
	- Cluster summary (clusters_summary.tsv)
	- Tab-delimited cluster summary table (clusters_summary_table.tsv)
	- Subdirectory per specified gap allowance (gap_#) containing:
		- CLUSTERS:
			- Round-robin reconstructed syntenic clusters for each species
		- PAIRS:
			- Round-robin identified gene pairs for each species

### <b>Step by step examples</b>

#### Example 1 - <i>Cryptococcus</i>
Below is a quick example describing how to compare two genomes from [<i>Cryptococcus neoformans</i> var. <i>neoformans</i> JEC21](https://pubmed.ncbi.nlm.nih.gov/15653466/) and [<i>Cryptococcus gattii</i> WM276](https://pubmed.ncbi.nlm.nih.gov/21304167/) using annotation data available in public databases.

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

##### Running SYNY and plotting comparisons with Circos:
```Bash
SYNY=~/SYNY_RESULTS      ## Replace by desired SYNY output directory

run_syny.pl \
  -a $DATA/*.gbff.gz \
  -g 0 1 5 \
  -e 1e-10 \
  -r JEC21 \
  -o $SYNY \
  --circos \
  --circos_prefix WM276_vs_JEC21
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


##### Example of an image generated with Circos and SYNY from shared proteins clusters:
<p align="left">
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/WM276_vs_JEC21.syny.normal.png">
</p>

In this figure, nucleotides biases are plotted in the concentric rings (from outer to inner rings):
- AT and GC nucleotide biases (grey and red lines)
- GT and AC nucleotide biases (blue and green lines)
- GA and CT nucleotide biases (purple and yellow lines)

Syntenic blocks identified by SYNY are indicated by ribbons. These ribbons are color-coded based on the chromosomes/contigs present in the reference genome used. The reference genome can be specified with the `--ref` commmand line switch. If omitted, the  first genome encountered alphabetically will be used as the default reference.

Data and configuration files for these plots are located in the CIRCOS/ subdirectory.

##### Example of an image generated with Circos and SYNY from shared proteins clusters (using the inverted karyotype):
<p align="left">
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/WM276_vs_JEC21.syny.inverted.png">
</p>

In the inverted karyotype image, the order of the karyotype(s) to be compared to the reference one are reversed. This option can be useful when comparing genomes whose chromosomes have been assigned similar numbers based on various inference methods (this does not appear to be the case in the above example). In such instances, inverting the karyotypes can help improve figure legibility.

##### Example of an image generated with Circos and SYNY from pairwise genome alignments:
<p align="left">
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/WM276_vs_JEC21.paf.normal.png">
</p>

In pairwise genome alignments, repetitive regions (such as telomeres/subtelomeres) can produce more than one alignment for a given locus. In the above figure, a bit of extra noise is added to the figure (as thin criss-crossing lines) due to these repetitive segments. As a rule of thumb, repetitive segments are easier to spot in dotplot-like figures (see dotplot section below).

##### Example of a barplot generated from the minimap2 pairwize alignments (PAF) files (using defaults settings):
<p align="left">
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/WM276_vs_JEC21.barplot.Spectral.png">
</p>

If pairwise genome alignments are performed with [minimap2](https://github.com/lh3/minimap2), barplots will be generated from the minimap2-generated PAF alignment files with [paf_to_barplot.py](https://github.com/PombertLab/SYNY/blob/main/paf_to_barplot.py) and [matplotlib](https://matplotlib.org/).

In these plots, colinear regions found between the compared genomes are highlighted by colored rectangles. By default, these rectangles are color-coded based on the contigs/chromosomes of the query. The above barplot image was generated using the Spectral color palette from [seaborn](https://seaborn.pydata.org/tutorial/color_palettes.html) (set as default in SYNY). This palette can be replaced using the `--palette` command line switch; <i>e.g.</i> `--palette husl` (see this [URL](https://www.practicalpythonfordatascience.com/ap_seaborn_palette) for a detailed list of available palettes).

If desired, the barplots can instead be generated using a single monochromatic color with the `--monobar` command line switch; <i>e.g.</i> `--monobar red`.

By default, the barplots are formatted for a widescreen (landscape) output (width/height ratio: 19.2/10.8). This ratio can be adjusted with the `--height` and `--width` command line switches.

##### Example of a dotplot generated from the minimap2 pairwize alignments (PAF) files (using defaults settings):
<p align="left">
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/WM276_vs_JEC21.1e5.blue.png">
</p>

If pairwise genome alignments are performed with [minimap2](https://github.com/lh3/minimap2), by default, dotplot-like scatter plots will be generated from the minimap2-generated PAF alignment files with [paf_to_dotplot.py](https://github.com/PombertLab/SYNY/blob/main/paf_to_dotplot.py) and [matplotlib](https://matplotlib.org/). If desired, this step can be skipped entirely with the `--no_dotplot` command line switch.

In these plots, each chromosome/contig from the query is plotted as a column (x-axis) against each chromosome/contig from the subject (y-axis). In the above example, a total of 196 subplots (14 x 14 chromosomes) are plotted using matplotlib's [subplot](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.subplots.html) function. In these plots, matches identified with minimap2 are scatter-plotted using the selected color (defaut: blue), with repeated loci indicated by the presence of matches across two or more contigs/chromosomes.

Because large numbers quickly overlap in the small subplots, to improve legibility, numbers in the x and y axes are scaled down using a desired scale with the `--multi` command line switch (default: 1e5). For example:
```bash
1e3: n x 1000 bp
1e4: n x 10000 bp
1e5: n x 100000 bp
1e6: n x 1000000 bp
```

Alternatively, ticks/numbers in the x and y axes can be turned off with the `--noticks` command line option.

By default, the dotplots are formatted for a widescreen (landscape) output (width/height ratio: 19.2/10.8). This ratio can be adjusted with the `--height` and `--width` command line switches (shared with [paf_to_barplot.py](https://github.com/PombertLab/SYNY/blob/main/paf_to_barplot.py)).

The default monochromatic color (blue) can be changed with the `--color` option (<i>e.g.</i> `--color red`). If desired, dotplots can instead be color-coded based on the query contigs/chromosomes with the the `--dotpalette` option; e.g. `--dotpalette inferno` (see this [URL](https://www.practicalpythonfordatascience.com/ap_seaborn_palette) for a detailed list of available palettes).

Note that plotting large genomes can quickly eat up a lot of memory. When running out of memory, the process will be terminated automatically before a PNG image can be produced.

#### Example 2 - <i>Encephalitozoon</i>
Below is a quick example describing how to compare a total of three telomere-to-telomere (T2T) genomes from <i>Encephalitozoon</i> species [<i>E. intestinalis</i> ATCC 50506](https://pubmed.ncbi.nlm.nih.gov/37142951/), [<i>E. hellem</i> ATCC 50604](https://pubmed.ncbi.nlm.nih.gov/37142951/), and [<i>E. cuniculi</i> ATCC 50602](https://pubmed.ncbi.nlm.nih.gov/37142951/) using annotation data available in public databases.

##### Downloading annotation data from GenBank (NCBI):
```bash
DATA=~/ENCE
mkdir -p $DATA

## Encephalitozoon intestinalis ATCC 50506 telomere-to-telomere (T2T) genome
curl \
  -L https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/024/399/295/GCA_024399295.1_ASM2439929v1/GCA_024399295.1_ASM2439929v1_genomic.gbff.gz \
  -o $DATA/intestinalis_50506.gbff.gz

## Encephalitozoon hellem ATCC 50604 telomere-to-telomere (T2T) genome
curl \
  -L https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/024/399/255/GCA_024399255.1_ASM2439925v1/GCA_024399255.1_ASM2439925v1_genomic.gbff.gz \
  -o $DATA/hellem_50604.gbff.gz

## Encephalitozoon cuniculi ATCC 50602 telomere-to-telomere (T2T) genome
curl \
  -L https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/027/571/585/GCA_027571585.1_ASM2757158v1/GCA_027571585.1_ASM2757158v1_genomic.gbff.gz \
  -o $DATA/cuniculi_50602.gbff.gz

```

##### Running SYNY and plotting comparisons with Circos:
```Bash
SYNY=~/SYNY_ENCE      ## Replace by desired SYNY output directory

run_syny.pl \
  -a $DATA/*.gbff.gz \
  -g 0 1 5 \
  -e 1e-10 \
  -r intestinalis_50506 \
  -o $SYNY \
  --circos \
  --circos_prefix encephalitozoon
```

##### Example of clusters identified with SYNY
```Bash
head -n 32 $SYNY/SYNTENY/clusters_summary.tsv

##### cuniculi_50602_vs_hellem_50604; Gap = 0 #####
  Total number of proteins in clusters: 1827
  # of clusters:        86
  Longest:      109
  Shortest:     2
  Average cluster size: 21
  Median cluster size:  12
  N50:  38
  N75:  23
  N90:  10

##### cuniculi_50602_vs_hellem_50604; Gap = 1 #####
  Total number of proteins in clusters: 1834
  # of clusters:        24
  Longest:      184
  Shortest:     2
  Average cluster size: 76
  Median cluster size:  48
  N50:  148
  N75:  128
  N90:  48

##### cuniculi_50602_vs_hellem_50604; Gap = 5 #####
  Total number of proteins in clusters: 1838
  # of clusters:        19
  Longest:      194
  Shortest:     2
  Average cluster size: 97
  Median cluster size:  92
  N50:  160
  N75:  146
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

#### Example 3 - <i>Encephalitozoon</i> with custom colors
Custom colors for [Circos](https://circos.ca/) plots can be loaded directly from tab-delimited text files containing color names and their associated RGB values (see [custom_color_1.txt](https://github.com/PombertLab/SYNY/blob/main/Circos/custom_color_1.txt) for an example). A few custom presets are also available to use in SYNY.

##### Running SYNY with custom colors loaded from a tab-delimited file:
```Bash
SYNY=~/SYNY_ENCE_CC         ## Replace by desired SYNY output directory
COLORS=~/custom_color_2.txt ## Replace by desired custom color file

run_syny.pl \
  -a $DATA/*.gbff.gz \
  -g 0 1 5 \
  -e 1e-10 \
  -r intestinalis_50506 \
  -o $SYNY \
  --custom_file $COLORS \
  --circos \
  --circos_prefix encephalitozoon_cc
```

##### Example of an image generated with Circos and SYNY comparing a total of 3 genomes (using a custom color set):
<p align="left">
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/encephalitozoon_cc.png">
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
SYNY=~/SYNY_ENCE_CC_PRESET   ## Replace by desired SYNY output directory

run_syny.pl \
  -a $DATA/*.gbff.gz \
  -g 0 1 5 \
  -e 1e-10 \
  -r intestinalis_50506 \
  -o $SYNY \
  --circos \
  --circos_prefix encephalitozoon_blues \
  --custom_preset blues
```

##### Example of an image generated with the blues custom color preset:
<p align="left">
  <img src="https://github.com/PombertLab/SYNY/blob/main/Images/encephalitozoon_blues.png">
</p>


## <b>References</b>
[Sensitive protein alignments at tree-of-life scale using DIAMOND](https://www.nature.com/articles/s41592-021-01101-x). Buchfink B, Reuter K, Drost HG. <b>Nature Methods.</b> 18, 366–368 (2021). doi:10.1038/s41592-021-01101-x

[Minimap2: pairwise alignment for nucleotide sequences](https://pubmed.ncbi.nlm.nih.gov/29750242/). Li H. <b>Bioinformatics.</b> 2018 Sep 15;34(18):3094-3100. doi: 10.1093/bioinformatics/bty191.

[Circos: an information aesthetic for comparative genomics](https://pubmed.ncbi.nlm.nih.gov/19541911/). Krzywinski M, Schein J, Birol I, Connors J, Gascoyne R, Horsman D, Jones SJ, Marra MA. <b>Genome Res.</b> 2009 Sep;19(9):1639-45. doi:10.1101/gr.092759.109

[The genome of the basidiomycetous yeast and human pathogen <i>Cryptococcus neoformans</i>](https://pubmed.ncbi.nlm.nih.gov/15653466/). Loftus BJ, Fung E, Roncaglia P, Rowley D, Amedeo P, Bruno D, Vamathevan J, Miranda M, Anderson IJ, Fraser JA, Allen JE, Bosdet IE, Brent MR, Chiu R, Doering TL, Donlin MJ, D'Souza CA, Fox DS, Grinberg V, Fu J, Fukushima M, Haas BJ, Huang JC, Janbon G, Jones SJ, Koo HL, Krzywinski MI, Kwon-Chung JK, Lengeler KB, Maiti R, Marra MA, Marra RE, Mathewson CA, Mitchell TG, Pertea M, Riggs FR, Salzberg SL, Schein JE, Shvartsbeyn A, Shin H, Shumway M, Specht CA, Suh BB, Tenney A, Utterback TR, Wickes BL, Wortman JR, Wye NH, Kronstad JW, Lodge JK, Heitman J, Davis RW, Fraser CM, Hyman RW. <b>Science.</b> 2005 Feb 25;307(5713):1321-4. doi: 10.1126/science.1103773. Epub 2005 Jan 13. PMID: 15653466; PMCID: PMC3520129.

[Genome variation in <i>Cryptococcus gattii</i>, an emerging pathogen of immunocompetent hosts](https://pubmed.ncbi.nlm.nih.gov/21304167/). D'Souza CA, Kronstad JW, Taylor G, Warren R, Yuen M, Hu G, Jung WH, Sham A, Kidd SE, Tangen K, Lee N, Zeilmaker T, Sawkins J, McVicker G, Shah S, Gnerre S, Griggs A, Zeng Q, Bartlett K, Li W, Wang X, Heitman J, Stajich JE, Fraser JA, Meyer W, Carter D, Schein J, Krzywinski M, Kwon-Chung KJ, Varma A, Wang J, Brunham R, Fyfe M, Ouellette BF, Siddiqui A, Marra M, Jones S, Holt R, Birren BW, Galagan JE, Cuomo CA. <b>mBio</b>. 2011 Feb 8;2(1):e00342-10. doi: 10.1128/mBio.00342-10. PMID: 21304167; PMCID: PMC3037005.

[Telomere-to-Telomere genome assemblies of human-infecting <i>Encephalitozoon</i> species](https://pubmed.ncbi.nlm.nih.gov/37142951/). Mascarenhas Dos Santos AC, Julian AT, Liang P, Juárez O, Pombert JF. <b>BMC Genomics</b>. 2023 May 4;24(1):237. doi: 10.1186/s12864-023-09331-3. PMID: 37142951; PMCID: PMC10158259.