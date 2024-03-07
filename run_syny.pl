#!/usr/bin/env perl
# Pombert lab, 2022

my $name = 'run_syny.pl';
my $version = '0.5.5';
my $updated = '2024-03-07';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;
use Cwd qw(abs_path);
use File::Path qw(make_path);

my $usage = <<"EXIT";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Runs all the components of the SYNY pipeline and
		creates configuration files for Circos plots

REQS		PerlIO::gzip
		Diamond - https://github.com/bbuchfink/diamond

USAGE		${name} \\
		  -a *.gbff \\
		  -e 1e-10 \\
		  -g 0 1 5 \\
		  -o SYNY \\
		  -r CCMP1205

OPTIONS (MAIN):
-a (--annot)	GenBank GBF/GBFF Annotation files (GZIP files are supported)
-e (--evalue)	BLAST evalue cutoff [Default = 1e-10]
-g (--gaps)	Allowable number of gaps between pairs [Default = 0]
-o (--outdir)	Output directory [Default = SYNY]

OPTIONS (PLOTS): ##### Requires Circos - http://circos.ca/ #####
-c (--circos)	Generate Circos plots automatically
-circos_prefix	Desired Circos plot prefix [Default: circos]
-r (--ref)	Genome to use as reference (defaults to first one alphabetically if none provided)
-u (--unit)	Size unit (Kb or Mb) [Default: Mb]
-custom_file	Load custom colors from file
-list_preset	List available custom color presets
-custom_preset	Use a custom color preset; e.g.
		# chloropicon - 20 colors - Lemieux et al. (2019) https://pubmed.ncbi.nlm.nih.gov/31492891/
		# encephalitozoon - 11 colors - Pombert et al. (2012) https://pubmed.ncbi.nlm.nih.gov/22802648/
EXIT

die ("\n$usage\n") unless (@ARGV);

my @commands = @ARGV;

my @annot_files;
my $evalue = '1e-10';
my @gaps;
my $outdir = 'SYNY';
my $reference;
my $unit = 'Mb';
my $circos;
my $circos_prefix = 'circos';
my $custom_file;
my $custom_colors;
my $list_preset;
my @formats;

GetOptions(
	'a|annot=s@{1,}' => \@annot_files,
	'e|evalue=s' => \$evalue,
	'g|gaps=s{0,}' => \@gaps,
	'o|outdir=s' => \$outdir,
	'r|ref|reference=s' => \$reference,
	'u|unit=s' => \$unit,
	'c|circos' => \$circos,
	'circos_prefix=s' => \$circos_prefix,
	'custom_file=s' => \$custom_file,
	'custom_preset=s' => \$custom_colors,
	'list_preset'	=> \$list_preset,
	'f|format=s{1,}' => \@formats,
);

unless(@gaps){
	@gaps = (0);
}

###################################################################################################
## Output directory creation and setup
###################################################################################################

my ($script,$path) = fileparse($0);
my $circos_path = "$path/Circos";

### List available custom presets for Circos, then exit
if ($list_preset){

	system("
		$circos_path/nucleotide_biases.pl \\
		--list_preset
	") == 0 or checksig();
	exit;

}

### Otherwize, create output dir

unless(-d $outdir){
	make_path($outdir,{mode => 0755}) or die "Can't create $outdir: $!\n";
}

## Convert $outdir to absolute path: required for Circos files
if ($outdir !~ /^\//){
	$outdir = abs_path($outdir);
}

my $list_dir = "$outdir/LISTS";
my $prot_dir = "$outdir/PROT_SEQ";
my $diamond_dir = "$outdir/DIAMOND";
my $db_dir = "$diamond_dir/DB";
my $shared_dir = "$outdir/SHARED";
my $synteny_dir = "$outdir/SYNTENY";
my $conserved_dir = "$outdir/CONSERVED";

my @outdirs = ($list_dir,$prot_dir,$diamond_dir,$db_dir,$synteny_dir);

foreach my $dir (@outdirs){
	unless (-d $dir){
		make_path($dir,{mode=>0755}) or die "Can't create $dir: $!\n";
	}
}

## Logs
my $log_file = $outdir.'/'.'syny.log';
my $log_err = $outdir.'/'.'error.log';

open LOG, '>', $log_file or die "Can't create $log_err: $!\n";
open ERROR, '>', $log_err or die "Can't create $log_err: $!\n";

my $time = localtime();
my $tstart = time();

print LOG "SYNY started on: ".$time."\n";
print LOG "COMMAND: $0 @commands\n";

###################################################################################################
## Run list_maker.pl
###################################################################################################

print ERROR "\n### list_maker.pl ###\n";

system("
	$path/list_maker.pl \\
	  --input @annot_files \\
	  --outdir $outdir \\
	  2>> $outdir/error.log
") == 0 or checksig();

###################################################################################################
## Run get_homology.pl
###################################################################################################

my %linked_files;
my @prot_files;

opendir (FAA, "$outdir/PROT_SEQ") or die "\n\n[ERROR]\tCan't open $outdir/PROT_SEQ: $!\n\n";

while (my $file = readdir(FAA)){
	if ($file =~ /\.faa$/){
		push (@prot_files, "$outdir/PROT_SEQ/$file");
	}
}

link_files();

print ERROR "\n### get_homology.pl ###\n";

system("
	$path/get_homology.pl \\
	--input @prot_files \\
	--evalue $evalue \\
	--outdir $diamond_dir \\
	--shared $shared_dir \\
	--list $list_dir \\
	2>> $outdir/error.log
") == 0 or checksig();

###################################################################################################
## Run get_synteny.pl
###################################################################################################

print ERROR "\n### get_synteny.pl ###\n";

my $clu_sum_file = $synteny_dir.'/'.'clusters_summary.tsv';
if (-f $clu_sum_file){
	system("rm $clu_sum_file") == 0 or checksig();
}

foreach my $annot_file_1 (sort(@annot_files)){

	my ($file_name_1,$dir,$ext) = fileparse($annot_file_1,'\..*');
	my $linked_file_1 = $linked_files{$annot_file_1};
	print "\nIdentifying synteny between $file_name_1\n"; 

	foreach my $annot_file_2 (sort(@annot_files)){
		if($annot_file_1 ne $annot_file_2){

			my ($file_name_2,$dir,$ext) = fileparse($annot_file_2,'\..*');
			my $linked_file_2 = $linked_files{$annot_file_2};
			print "\t$file_name_2\n";

			foreach my $gap (@gaps){
				system("
					$path/get_synteny.pl \\
					--query_list $list_dir/$file_name_1.list \\
					--query_blast $diamond_dir/${file_name_1}_vs_${file_name_2}.diamond.6 \\
					--subject_list $list_dir/$file_name_2.list \\
					--subject_blast $diamond_dir/${file_name_2}_vs_${file_name_1}.diamond.6 \\
					--gap $gap \\
					--outdir $synteny_dir/gap_$gap \\
					--sumdir $synteny_dir \\
					2>> $outdir/error.log
				") == 0 or checksig();
			}
		}
		else{
			print "\t. . .\n";
		}
	}

}

### Grabbing total number of entries in .list files
my %gene_count;
opendir (LISTS, "$outdir/LISTS") or die "\n\n[ERROR]\tCan't open $outdir/LISTS: $!\n\n";

while (my $file = readdir(LISTS)){

	if ($file =~ /\.list$/){

		my $list_file = $outdir.'/LISTS/'.$file;
		open TMP, '<', $list_file or die "$!\n";
		my ($basename,$path) = fileparse($file);
		$basename =~ s/\.list$//;

		my $counter = 0;
		while (my $line = <TMP>){
			chomp $line;
			if ($line ne ''){
				$counter++;
			}
		}
		$gene_count{$basename} = $counter;

	}
}

### Create a cluster summary table
my $clu_sum_table = $synteny_dir.'/'.'clusters_summary_table.tsv';
open CLU, '<', $clu_sum_file or die "Can't read $clu_sum_file: $!\n";
open TSV, '>', $clu_sum_table or die "Can't create $clu_sum_table: $!\n";

print TSV '### Query'."\t";
print TSV 'Total # proteins'."\t";
print TSV 'Allowed Gaps'."\t";
print TSV 'Total # proteins in clusters'."\t";
print TSV '% of proteins in clusters'."\t";
print TSV '# of clusters'."\t";
print TSV 'Longest'."\t".'Shortest'."\t";
print TSV 'Average'."\t".'Median'."\t";
print TSV 'N50'."\t".'N75'."\t".'N90'."\n";

my %cluster_metrics;
my $clu_query;
my $qgap;

while (my $line = <CLU>){
	chomp $line;
	if ($line =~ /##### (\w+); Gap = (\d+) #####/){
		$clu_query = $1;
		$qgap = $2;
	}
	elsif ($line =~ /Total.*\t(\S+)$/){
		$cluster_metrics{$clu_query}{$qgap}{'total'} = $1;
	}
	elsif ($line =~ /# of clusters.*\t(\S+)$/){
		$cluster_metrics{$clu_query}{$qgap}{'total_clu'} = $1;
	}
	elsif ($line =~ /Longest.*\t(\S+)$/){
		$cluster_metrics{$clu_query}{$qgap}{'longest'} = $1;
	}
	elsif ($line =~ /Shortest.*\t(\S+)$/){
		$cluster_metrics{$clu_query}{$qgap}{'shortest'} = $1;
	}
	elsif ($line =~ /Average.*\t(\S+)$/){
		$cluster_metrics{$clu_query}{$qgap}{'average'} = $1;
	}
	elsif ($line =~ /Median.*\t(\S+)$/){
		$cluster_metrics{$clu_query}{$qgap}{'median'} = $1;
	}
	elsif ($line =~ /N50.*\t(\S+)$/){
		$cluster_metrics{$clu_query}{$qgap}{'n50'} = $1;
	}
	elsif ($line =~ /N75.*\t(\S+)$/){
		$cluster_metrics{$clu_query}{$qgap}{'n75'} = $1;
	}
	elsif ($line =~ /N90.*\t(\S+)$/){
		$cluster_metrics{$clu_query}{$qgap}{'n90'} = $1;
	}
}

foreach my $query (sort (keys %cluster_metrics)){
	foreach my $gap (sort {$a <=> $b}(keys %{$cluster_metrics{$query}})){

		print TSV $query."\t";

		my ($ref) = $query =~ /^(\w+)\_vs/;
		print TSV $gene_count{$ref}."\t";

		print TSV $gap."\t";

		print TSV $cluster_metrics{$query}{$gap}{'total'}."\t";
		my $percent = ($cluster_metrics{$query}{$gap}{'total'}/$gene_count{$ref})*100;
		$percent = sprintf("%.2f", $percent);
		print TSV $percent."\t";

		print TSV $cluster_metrics{$query}{$gap}{'total_clu'}."\t";
		print TSV $cluster_metrics{$query}{$gap}{'longest'}."\t";
		print TSV $cluster_metrics{$query}{$gap}{'shortest'}."\t";
		print TSV $cluster_metrics{$query}{$gap}{'average'}."\t";
		print TSV $cluster_metrics{$query}{$gap}{'median'}."\t";
		print TSV $cluster_metrics{$query}{$gap}{'n50'}."\t";
		print TSV $cluster_metrics{$query}{$gap}{'n75'}."\t";
		print TSV $cluster_metrics{$query}{$gap}{'n90'}."\n";

	}
}

close CLU;
close TSV;

###################################################################################################
## Run id_conserved_regions.pl
###################################################################################################

print ERROR "\n### id_conserved_regions.pl ###\n";

system("
	$path/id_conserved_regions.pl \\
	--lists $list_dir \\
	--blasts $diamond_dir \\
	--outdir $conserved_dir \\
	2>> $outdir/error.log
") == 0 or checksig();

###################################################################################################
## Create links files, karyotypes and nucleotide biases for Circos
###################################################################################################

## Circos links files
print ERROR "\n### clusters2links.pl ###\n";

my $cluster_dir = "$outdir/SYNTENY";

my @cluster_dirs;
opendir (CDIR, $cluster_dir) or die "Can't open $cluster_dir: $!\n";
while (my $dname = readdir(CDIR)) {
	if (-d "$cluster_dir/$dname"){
		unless (($dname eq '.') or ($dname eq '..')){
			push (@cluster_dirs, "$cluster_dir/$dname");
		}
	}
}

foreach my $cluster_subdir (@cluster_dirs){
	my ($basedir) = fileparse($cluster_subdir);
	system("
		$circos_path/clusters2links.pl \\
		--cluster $cluster_subdir/CLUSTERS/*.clusters \\
		--list $list_dir/*.list \\
		--outdir $outdir/CIRCOS \\
		2>> $outdir/error.log
	") == 0 or checksig();
}

## Karyotypes and nucleotide biases
print ERROR "\n### nucleotide_biases.pl ###\n";

# Option flags
my $ref = '';
if ($reference){
	$ref = "--reference $reference";
}

my $custom_cc_file = '';
if ($custom_file){
	$custom_cc_file = "--custom_file $custom_file";
}

my $custom_cc = '';
if ($custom_colors){
	$custom_cc = "--custom_preset $custom_colors";
}

my $unit_size = '';
if ($unit){
	$unit_size = "--unit $unit";
}

my $gap = $gaps[0];

# Running nucleotide_biases.pl
system("
	$circos_path/nucleotide_biases.pl \\
	--outdir $outdir/CIRCOS \\
	--fasta $outdir/GENOME/*.fasta \\
	--winsize 10000 \\
	--step 5000 \\
	--gap $gap \\
	$ref \\
	$unit_size \\
	$custom_cc_file \\
	$custom_cc \\
	2>> $outdir/error.log
") == 0 or checksig();

# Running Circos
if ($circos){

	my $normal_circos = "$circos_prefix.normal.png";
	my $invert_circos = "$circos_prefix.inverted.png";

	print "\nRunning Circos on genotype (normal): $normal_circos\n\n";

	system ("
	  circos \\
	  -conf $outdir/CIRCOS/concatenated/concatenated.conf \\
	  -outputdir $outdir \\
	  -outputfile $normal_circos");

	print "\nRunning Circos on genotype (inverted): $invert_circos\n\n";

	system ("
	  circos \\
	  -conf $outdir/CIRCOS/concatenated/concatenated.inverted.conf \\
	  -outputdir $outdir \\
	  -outputfile $invert_circos");
	
}

###################################################################################################
## Completion
###################################################################################################

my $end_time = localtime();
my $run_time = time - $tstart;

print LOG "SYNY completed on: ".$end_time."\n";
print LOG "Runtime: ".$run_time." seconds\n";

###################################################################################################
## Subroutines
###################################################################################################

sub link_files {

	if(scalar(@annot_files) != scalar(@prot_files)){
		die("The number of annotation files (".scalar(@annot_files).") does not equal the number of protein files (".scalar(@prot_files).")\n");
	}

	foreach my $annot_file (@annot_files){
		my ($annot_name,$dir,$ext) = fileparse($annot_file,"\..*");
		INNER: foreach my $prot_file (@prot_files){
			my ($prot_name,$dir,$ext) = fileparse($prot_file,"\..*");
			if ($prot_name eq $annot_name){
				$linked_files{$annot_file} = $prot_file;
				last INNER;
			}
		}
	}

	foreach my $annot_file (@annot_files){
		unless($linked_files{$annot_file}){
			print("[W]  No matching protein file was found for $annot_file. It will be skipped as a consequence.\n");
		}
	}

}

sub checksig {

	my $exit_code = $?;
	my $modulo = $exit_code % 255;

	if ($modulo == 2) {
		print "\nSIGINT detected: Ctrl+C => exiting...\n\n";
		exit(2);
	}
	elsif ($modulo == 131) {
		print "\nSIGTERM detected: Ctrl+\\ => exiting...\n\n";
		exit(131);
	}

}