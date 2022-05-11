The SYNY pipeline investigates synteny between species by reconstructing protein clusters from gene pairs generated from [DIAMOND](https://github.com/bbuchfink/diamond) homology searches.

## <b>Introduction</b>
#### <b>What is synteny?</b>
Synteny is the measure of similarity in gene organization between two species. Through time, genomes can be reorganized by large-scale events such as recombination, but also smaller events like -- but not exhaustively -- duplication, translocation, deletion, or inversion.
#### <b>Why use synteny?</b>
Synteny can be used in a variety of ways, making it a multi-tool analysis. Some examples of usage:<br>
- Evolutionary distance analysis
	- Species sharing similar/identical reorganization events are likley more closely related than species that do not share these events.
- Genome annotation correction
	- Identify mispredicted genes

## <b>Requirements</b>
- [Perl5](https://www.perl.org/)
- [Python-3.7+](https://www.python.org/downloads/release/python-370/)
	- [matplotlib](https://matplotlib.org/)

## <b>Using SYNY</b>
The full SYNY pipeline can be run with [run_syny.pl](), that will:
1. Extract protein sequences from provided annotation files (currently supported: .gbf, .gbff, .prot, .embl).
2. Preform round-robin [DIAMOND](https://github.com/bbuchfink/diamond) BLASTP homology searches.
3. Generate gene pairs and reconstruct gene clusters.
4. Identify sample-wide conserved genes, as well as missing or unique proteins for each species.

SYNY can be run utilizing the master script as follows:<br>
```
run_SYNY.pl \
  -a *.gbff \
  -g 5 \
  -o Chloropicon_synteny
```
Options for run_SYNY.pl are:
```
-a (--annot)	Annotation files (Supported files: gbff, gff, embl)
-e (--evalue)	BLAST evalue cutoff [Default = 1e-10]
-g (--gaps)		Allowable number of gaps between pairs [Default = 0]
-o (--outdir)	Output directory [Default = SYNY]
```
The output directory will be structured as follows: 
```
drwxr-xr-x. 2 julian julian 4.0K May 10 18:18 CONSERVED
drwxr-xr-x. 3 julian julian 4.0K May 10 18:18 DIAMOND
drwxr-xr-x. 3 julian julian 4.0K May 10 19:29 LISTS
drwxr-xr-x. 2 julian julian 4.0K May 10 19:29 PROT_SEQ
drwxr-xr-x. 3 julian julian 4.0K May 10 18:18 SYNTENY
```
- CONSERVED:
	- List of proteins that are conserved within the sample, for each species (.conserved)
	- Summary of all homology results and conservation within the sample, for each species (.conserved_summary)
	- List of unique proteins for each species (.unique)
- DIAMOND:
	- DB
		- BLASTP databases
	- Round-robin BLASTP results (.diamond.6)
- LISTS:
	- Lists connecting accession numbers to locus tags (ALIASES)
	- Lists of protein coding genes with location details (.list)
- PROT_SEQ:
	- Protein sequences for each species (.faa)
- SYNTENY:
	- Directory per specified gap allowance (gap_#)
		- CLUSTERS:
			- Round-robin reconstructed syntetic clusters for each species
		- PAIRS:
			- Round-robin identified gene pairs for each species

## <b>SYNY step-by-step [Coming soon]</b>

## <b>References</b>
[Sensitive protein alignments at tree-of-life scale using DIAMOND](https://www.nature.com/articles/s41592-021-01101-x). <b>Buchfink B, Reuter K, Drost HG.</b> Nature Methods 18, 366–368 (2021). doi:10.1038/s41592-021-01101-x
