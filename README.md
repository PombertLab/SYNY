## <b>Synopsis</b>

The SYNY pipeline investigates gene colinearity (synteny) between genomes by reconstructing clusters from conserved pairs of protein-coding genes identified from [DIAMOND](https://github.com/bbuchfink/diamond) homology searches.

## <b>Introduction</b>
#### <b>What is synteny?</b>
Synteny is the measure of similarity in gene organization between two species. Through time, genomes can be reorganized by large-scale events such as recombination, but also smaller events like -- but not exhaustively -- duplication, translocation, deletion, or inversion.

#### <b>Why use synteny?</b>
Synteny inferences can be used to:
- Improve genome annotation (by helping to identify mispredicted and/or missing genes)
- Perform evolutionary distance analyses (<i>e.g.</i> species sharing similar/identical genome reorganization events are likely more closely related than species that do not share these events).

## <b>Requirements</b>
- [Diamond](https://github.com/bbuchfink/diamond) ## To perform homology searches
- [Perl5](https://www.perl.org/)
- [PerlIO::gzip](https://metacpan.org/pod/PerlIO::gzip)

##### <b>Optional</b>
- [Circos](https://circos.ca/) ## To generate circos plots

## <b>Using SYNY</b>
The SYNY pipeline can be run with [run_syny.pl](https://github.com/PombertLab/SYNY/blob/main/run_syny.pl), a master script that:
1. Extracts protein sequences from provided annotation files (currently supported: .gbf, .gbff, .prot, .embl).
2. Performs round-robin [DIAMOND](https://github.com/bbuchfink/diamond) BLASTP homology searches.
3. Generates gene pairs and reconstruct gene clusters.
4. Identifies sample-wide conserved genes, as well as missing or unique proteins for each species.
5. Generates configuration files/templates for plotting with [Circos](https://circos.ca/).

SYNY can be run utilizing the master script as follows:<br>
```
run_SYNY.pl \
  -a *.gbff \
  -g 5 \
  -o Chloropicon_synteny
```
Options for run_SYNY.pl are:
```
OPTIONS (MAIN):
-a (--annot)	Annotation files (Supported files: gbff, gff, embl)
-e (--evalue)	BLAST evalue cutoff [Default = 1e-10]
-g (--gaps)		Allowable number of gaps between pairs [Default = 0]
-o (--outdir)	Output directory [Default = SYNY]
-p (--prot)		Protein files # Generated automatically if GenBank files

OPTIONS (PLOTS): ##### Requires Circos - http://circos.ca/ #####
-r (--ref)		Genome to use as reference (defaults to first one
				alphabetically if none provided)
-u (--unit)		Size unit (Kb or Mb) [Default: Mb]
-c (--circos)	Generate Circos plots; currently buggy
				# works if run independently on configuration files
				# generated, e.g.: circos --conf concatenated.conf
-custom			Use custom color palette (20 colors) from Lemieux et al.:
				# https://pubmed.ncbi.nlm.nih.gov/31492891/
```
The output directory will be structured as follows: 
```
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
	- Contains tab-delimited lists of locus_tags and their products
- CIRCOS:
	- Contains configuration files for plotting with Circos
- CONSERVED:
	- List of proteins that are conserved within the sample, for each species (.conserved)
	- Summary of all homology results and conservation within the sample, for each species (.conserved_summary)
	- List of unique proteins for each species (.unique)
- DIAMOND:
	- DB
		- BLASTP databases
	- Round-robin BLASTP results (.diamond.6)
- GENOME:
	- Contains FASTA files from the genomes being investigated
- LISTS:
	- Lists connecting accession numbers to locus tags (ALIASES)
	- Lists of protein coding genes with location details (.list)
- PROT_SEQ:
	- Protein sequences for each species (.faa)
- SHARED:
	- Contains lists of all proteins and their top homologs (if any) in other species (.shared.tsv)
	- Contains lists of proteins that are unique to each species (.uniques.tsv)
- SYNTENY:
	- Directory per specified gap allowance (gap_#)
		- CLUSTERS:
			- Round-robin reconstructed syntetic clusters for each species
		- PAIRS:
			- Round-robin identified gene pairs for each species

## <b>SYNY step-by-step [Coming soon]</b>




## <b>References</b>
[Sensitive protein alignments at tree-of-life scale using DIAMOND](https://www.nature.com/articles/s41592-021-01101-x). <b>Buchfink B, Reuter K, Drost HG.</b> Nature Methods 18, 366–368 (2021). doi:10.1038/s41592-021-01101-x
