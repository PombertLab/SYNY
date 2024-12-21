# SYNY CHANGE LOG

## SYNY-v1.2a
- `get_paf.pl` now generates VCF files from minimap2 genome alignment (min. alignment lenght = 1000 bp)
- `nucleotide_biases.pl` now calculates GC and AT skews. Corresponding data files are located in the `PLOTS/CIRCOS_DATA/` subdirectory.
- GC/AT skews are now plotted automatically with Circos. If desired, these subplots can be turned off independently with the `--no_skews` option, or together with all nucleotide biases subplots (with `--no_ntbiases`).
- Added a simple Fasta + GFF3 to GBFF converter (`gff3_to_gbff.pl`) in the `Utils/` subdirectory. This tool was tested on NCBI GFF3 files and expects the GFF3 file(s) to include gene/mRNA/exon/CDS entries in the `type` column and the `ID` and `Parent` tags in the attributes column. It also expects the corresponding Fasta and GFF3 files to share the same prefixes (e.g. genome_1.fasta / genome_1.gff). The GBFF files thus created were designed to work with SYNY but do not adhere exactly to the GBFF format and may not work for other purposes.
- `list_maker.pl`/ `run_syny.pl`: GenBank Flat file format extensions (gbk, gb, gbf) are now recognized/accepted
- `check_mp_colors.py`: removed obsolete references to pylab
- Added `orient_fastas_to_reference.py` to the `Utils/` subdirectory. This script reorients contigs in FASTA file(s) based on BLASTN homology searches against a reference. This can be useful when working with newly assembled genomes.

## SYNY-v1.2
### Bugfix release

- Fixed concatenation issue with isoforms in `list_maker.pl`
- Fixed subranges issues in `list_maker.pl`
- Adjusted linearmap alpha value and edge color for readability in `linear_maps.py`
- Slightly reduced memory usage with matplotlib

## SYNY-v1.1b
### Bugfix release

- Fixed extra length issues with barplots, dotplots and linemaps. Code was missing a line.strip(). Issue created visual artefacts on barplots (longer frames).

##### Code cleanup:
- Added `--version` option for all scripts.
- Minor code cleanup / standardisation across scripts

## SYNY-v1.1a

##### Additions
- Added the `--include` option to select contigs by name from text file(s); one name per line
- Added the `--ranges` option to select contig subranges from text file(s); name start end
- Added the `--bpmode` option to generate pairwise (pair) and/or concatenated (cat) barplots. Possible values are `pair` (default), `cat`, and `all` (for both).
- Added the `--bclusters` option to color clusters by alternating colors in the barplots. The colors are not related within or between contigs, they are just used to highlight collinear chunks.
- Created `check_versions.pl` to summarize script versions; this information can now be displayed with `run_syny.pl --version`.

##### Bugfixes
- `list_maker.pl` now grabs GeneID tags if locus tags are absent from GBFF annotation files.
- Fixed .txt file extension + added a file size check to `paf_metrics.py`. Now skips plotting if file is empty.
- Fixed div by zero issue in `nucleotide_biases.pl`.
- Added a check to detect if annotations parsed are blank. `run_syny.pl` no longer crashes if annotations are blank when running gene cluster inferences. If blank, it now now skips this section automatically.
- Fixed perl env shebangs causing issues with conda
- Fixed wrong exit codes with readmes

##### Readme / logs
- Added section about memory usage with genome alignments
- Added mashmap barplot examples in the <i>Encephalitozoon</i> section
- Added `changes.md` summarizing changes between versions
- Improved `syny.log` file.

## SYNY-v1.1

##### Additions:
- SYNY now generates linear maps (aka linemaps) from PAF files with `linear_maps.py`.
- Added support for MashMap3 genome alignments. Mashmap can be selected instead of minimap with `--aligner mashmap`. It runs in a smaller memory footprint than minimap (if using its default percentage identity of 85%). It does not product exact alignments however.
- Added the option to exclude contigs by name matching regular expression(s): e.g. `--exclude '^AUX' '^CPGT'`.
- Added an alternate SYNY installation method that does not require sudo privileges by leveraging conda packages.

##### Fixes:
- Fixed the `The number of annotation files (2) does not equal the number of protein files (1)` error => rewrote the corresponding segment and removed the obsoleted subroutine.
- Fixed the unreliable $diamond_check in `get_homology.pl` (i.e. replaced <i>which</i> by <i>command -v</i>).
- Changed Perl dependency Roman => Text::Roman in `nucleotide_biases.pl`.

## SYNY-v1.0b

- run_syny.pl options can now be set from a configuration file (requires Getopt::ArgvFile); e.g. `run_syny.pl @commands.conf`
- Added the Getopt::ArgvFile dependency to `setup_syny.pl` => `sudo cpanm Getopt::ArgvFile`
- Added a minimum contig size option + set defaults to all contigs, i.e. (`--minsize 1`)
- Added a matplotlib color palette check before computations so that plots won't crash if the color palette entered does not exist

## SYNY-v1.0a

- Added `--hfsize`, `--hmin`, `--hmax` and `--hauto` options to heatmaps
- Added more options to the Circos `--labels `command line switch. Possible values are now: `mixed`, `roman`, `arabic` and `names`
- Added `--pthreads` option to set the limit of plotting instances to run in parallel (in case each plot eats up too much RAM); defaults to the value set by `--threads` if omitted.
- Added SVG output to `paf_metrics.py`
- Set fonts as editable in SVG output files
- Removed unnecessary border frames from barplots
- Fixed ambiguous heatmap titles
- Added an example script (`Arabidopsis.sh`) in `Examples/` to download two Arabidopsis genomes (~ 100-150 Mbp each) for testing purposes

## SYNY-v1.0
### Stable release

#### Misc:
- Fixed output directory bug in `run_syny.pl` when using a deep tree
- Fixed abs_path() issue in `setup_syny.pl` that caused incomplete paths in the output configuration file
- Created `check_mp_colors.py` to list/plot color palettes available on the system (Fedora 40/Ubuntu 22.04 matplotlib palettes are not the same - 170 vs. 166) + added color palette plot
(Images/python_color_palettes.png)
- Fixed out of bounds barplot legends
- Added font size options -`-bfsize`/`--dfsize `options for barplots/dotplots

#### Circos:
- Contigs from the reference genome are now visually distinct and are labelled by roman numerals. Other contigs are labelled by arabic numerals.
- Added `--orientation` option (possible values: `normal`, `inverted`, `both`) + removed the now obsoleted `--no_invert`/`--no_normal` options
- Added `--no_cticks` option to disable ticks in Circos plots.
- added `--no_ntbiases` option to disable nucleotide bias subplots.
- Changed the default Circos plot mode to pairwise (`--circos pair`); concatenated plots can take a while to compute and are not always useful.
- Circos figures are now plotted in `--orientation normal` by default instead of both normal/inverted => less wasteful.
- Renamed the `.genotype` files generated by SYNY as `.karyotype` to match the nomenclature used by Circos

# OLDER VERSIONS

## SYNY-v0.9e
#### Cleanup release:

- Fixed a bug that crashed `nucleotide_biases.pl` when the reference entered was not found. Now uses the first sequence alphabetically if the ref entered is not found.
- Created `fasta_to_gbff.pl` to convert FASTA sequences to GBFF files (without annotations); useful to compare newly assembled genomes using pairwise alignments
- Added `Alignments`, `Clusters`, `Plots`, and `Utils` subdirs to the git repository and moved scripts/data accordingly
- Added shell scripts to download the example annotation data from NCBI
- Improved/cleaned up README

## SYNY-v0.9d
#### Cleanup release:

- Sanitized output directory:
    - Regrouped subdirs by analysis (`ALIGNMENTS/`, `CLUSTERS/`) and moved content accordingly
    - Created `PLOTS/` subdir and moved all plots therein
    - Renamed the CIRCOS data folder as `CIRCOS_DATA/` for greater clarity
    - Created `SEQUENCES/` subdir to store genome and protein fasta files
- Restructured/cleaned up `run_syny.pl`
- Improved the output log (`syny.log`)

## SYNY-v0.9c

- Simplified default help message
- Added `--help` option => displays all command lines options
- Added `--no_clus` option => turns off gene cluster inferences
- Added `--no_circos`, `--no_barplot` and `--no_heatmap` options => skips the correspoding plots
- Now generates both pairwise and concatenated Circos plots; `--circos all` is set as default

## SYNY-v0.9b
#### Bugfix release:
- Fixed strandedness in `clusters_to_paf.pl`, which caused collinearity SNAFUs in dotplots generated from the corresponding PAF files (e.g. `.gap_5.paf`).

## SYNY-v0.9a
- Parallelized several processes (list creation, PAF metrics, barplots/dotplots/heatmaps)
- Fixed an issue with out-of-bounds links in Circos plots due to 1-based vs. 0-based PAF columns
- Standardized shell output + added progress counters

## SYNY-v0.9
- `run_syny.pl`: Circos plotting is now multithreaded (one plot per thread)
- `run_syny.pl`: Now `generates concatenated and/or pairwise circos plots with --circos cat`, `--circos pair`, `--circos plot`
- `run_syny.pl`: Changed default threads value to 16 ## Previously defaulted to 8
- `nucleotide_biases.pl`: Now generates concatenated, pairwise and single circos configuration files

## SYNY-v0.8e
- `get_synteny.pl`: Fixed issue with out-of-order locus_tags when inferring synteny from protein clusters
- `list_maker.pl`: Fixed issue with multiple isoforms sharing the same locus_tag in GenBank files
- Now generates heatmaps summarizing percentages of colinear bases between genomes (with `paf_to_heatmap.py`)
- Fixed title in protein cluster heatmaps
- `--threads` option now applies to minimap2 alignments and diamond homology searches

## SYNY-v0.8d
- `run_syny.pl` now generates Circos plots for all requested gap values if `--circos` is invoked (+ moved plotting to subs).
- Moved Circos plots to `CIRCOS_PLOTS/` subdirectory
- Standardized Circos / barplot / dotplot file names using the `.mmap`/`.gap_0` affixes
- Renamed Circos configuration filenames in the` CIRCOS/` subdirectory for greater clarity

## SYNY-v0.8c
- `paf_to_barplot.py` / `paf_to_dotplot.py` now load queries/subjects from fasta files => otherwise some queries/subjects are missing from PAF files when no match is found.
- Fixed relative positions in PAF files generated with `clusters_to_paf.pl`
- Dotplots are now generated both from minimap2 pairwise genome alignments (`.mmap.`) and from protein clusters found with SYNY (e.g. `.gap_0.`)
- Added minimap2 `--threads` option to `run_syny.pl` + set default value to 8

## SYNY-v0.8b
- Fixed memory usage in `paf_to_dotplot.py` + added SVG output. Now runs much faster too.
- Fixed memory leak in `paf_to_barplot.py`
- `clusters_to_paf.pl` now generates PAF files from clusters identified by SYNY
- Barplots are now generated from these PAF files and are identified with the gap affix, e.g. `.gap_0`.
- Barplots generated from minimap2 alignments are now labelled with the `.mmap.` affix

## SYNY-v0.8a
- Moved PAF to Circos links conversion to `paf2links.pl` subscript
- Added `--clusters` option to color ribbons in Circos plots by clusters instead of by contigs ## Useful when comparing bacterial genomes

## SYNY-v0.8
- Added installation script `setup_syny.pl`. Tested on Fedora, Ubuntu, Debian, Kali and openSUSE Tumbleweed Linux distributions.
- Changed default heatmap palette from `crest` to `winter_r`; crest was missing from seaborn in some Linux distros...
- Fixed issue with clustered dendrograms; cm.figure.suptitle => cm.fig.suptitle; .figure.subtitle was not recognized in all distros...

## SYNY-v0.7f
- Now generates matrices sumarizing percentages of colinear protein-coding genes for each gap value investigated: e.g. `SYNTENY/gap_0/matrix_gap_0.tsv`
- Rewrote `protein_cluster_hm.py` to generate heatmaps from these matrices by leveraging pandas dataframes
- `protein_cluster_hm.py` now generates clustered dendrograms in addition to standard heatmaps
- Heatmaps and clustered dendrograms are now also generated in SVG format
- Added Circos `--labels`, `--label_size` and `--label_font` options; contigs can now be labelled by their names with `--labels names`

## SYNY-v0.7e
- Added heatmaps displaying the percentages of proteins found in clusters between each pair of genomes (summarized in `SYNTENY/clusters_summary_table.tsv`)
- Fixed dotplot issue with unidimensional arrays.
- Fixed dotplot issue when the total number of subplots is 1; now generates a single plot instead of a subplot.

## SYNY-v0.7d
- Now calculates pairwise genome alignment metrics from minimap2 PAF files, summarizes them as scatter plots, and stores the results in the `ALIGNMENTS/METRICS` subdirectory (see `paf_metrics.py`).
- Minor README restructuring to improve readability

## SYNY-v0.7c
- Added a `--resume` option to skip previously computed minimap2 alignments ## Useful when optimizing barplots/dotplots
- Added dimensions to barplot/dotplot output file names ## To prevent overwriting previous files when optimizing plots
- Barplot/dotplot height/width options are now independent: `--height`/`--dheight` and `--width`/`--dwidth`

## SYNY-v0.7b
- Added preset option for minimap2 (`--asm 5`, `--asm 10` or `--asm 20`), default = off
- Added option to skip dotplots (`--no_dotplot`)
- Added option to adjust gaps in dotplots (`--wdis` / `--hdis`)
- Added options to adjust Circos ticks, ideograms, links and points per track max values

## SYNY-v0.7a
- Reduced memory usage with `paf_to_dotplots.py`
- Added color scheme to barplot/dotplot output files
- Misc bug fixes

## SYNY-v0.7
- Now generates barplots showing colinear blocks between compared genomes using a per contig/chromosome color palette (default) or using a monochrome color instead (with `--monobar blue`).
- Added a color palette option to dotplots (e.g. `--dotpalette inferno`)
- Minor code fixes

## SYNY-v0.6b
- Cleaner dotplots
- Minor code fixes
- Better readme

## SYNY-v0.6a
- Now generates pairwise genome alignment dotplots (in PNG format) from the minimap2 PAF files (using mathplotlib)
- Minor code fixes

## SYNY-v0.6
- Added minimap2 dependency: now generates pairwise genome alignments in MAF, PAF and ALN (BLAST-like) formats
- Now generates colinearity plots from pairwise genome alignments (PAF) as well as from conserved protein clusters (SYNY)

## SYNY-v0.5c
- Circos plots are now generated properly

## SYNY-v0.5b
- simplified output: merged annotation/feature lists
- removed obsoleted options; support for GFF/EMBL files was broken

## SYNY-v0.5a
- Fixed crash on Ubuntu 22.04 with `get_homology.pl`
- Added Circos installation HOWTO from its tarball archive

##### Note: using `apt install circos` to install Circos on Ubuntu does not install its configuration files in the proper relative paths.

## SYNY-v0.5
- Initial release with code cleaned up a bit. Should work as intended.