#!/usr/bin/perl
# Pombert lab, 2022

my $name = 'run_syny.pl';
my $version = '0.5.2';
my $updated = '2023-09-25';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;
use File::Path qw(make_path);

my $usage = <<"EXIT";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Runs all the components of the SYNY pipeline step-by-step.

REQS		PerlIO::gzip

USAGE		${name} \\
		  -a *.gbff \\
		  -g 5 \\
		  -e 1e-10 \\
		  -r CCMP1205 \\
		  -o Chloropicon_synteny

OPTION
-a (--annot)	Annotation files (Supported files: gbff, gff, embl)
-e (--evalue)	BLAST evalue cutoff [Default = 1e-10]
-g (--gaps)	Allowable number of gaps between pairs [Default = 0]
-o (--outdir)	Output directory [Default = SYNY]
-r (--reference)	Genome to use as reference [for Circos]
-c (--circos)	Plot Circos images
-p (--prot)	Protein files ## Generated automatically from gbff files
EXIT

## No yet implemented:
# -f (--format)	Figure format(s) [Default: .svg]
# .png,.eps,.jpg,.jpeg,.pdf,.pgf,.ps,.raw,.rgba,.svg,.svgz,.tif,.tiff
# -n (--name)	Figure file prefix [Default: PHYLOGENETIC_FIGURE]

die ("\n$usage\n") unless (@ARGV);

my @annot_files;
my @prot_files;
my $evalue = '1e-10';
my @gaps;
my $outdir = 'SYNY';
my $reference;
my $circos;
my @formats;

GetOptions(
	'a|annot=s@{1,}' => \@annot_files,
	'p|proteins=s@{1,}' => \@prot_files,
	'e|evalue=s' => \$evalue,
	'g|gaps=s{0,}' => \@gaps,
	'o|outdir=s' => \$outdir,
	'r|reference=s' => \$reference,
	'c|circos' => \$circos,
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

unless(-d $outdir){
	make_path($outdir,{mode => 0755}) or die "Can't create $outdir: $!\n";
}

my $list_dir = "$outdir/LISTS";
my $prot_dir = "$outdir/PROT_SEQ";
my $annot_dir = "$outdir/ANNOTATIONS";
my $diamond_dir = "$outdir/DIAMOND";
my $db_dir = "$diamond_dir/DB";
my $synteny_dir = "$outdir/SYNTENY";
my $conserved_dir = "$outdir/CONSERVED";

my @outdirs = ($list_dir,$prot_dir,$annot_dir,$diamond_dir,$db_dir,$synteny_dir);

foreach my $dir (@outdirs){
	unless (-d $dir){
		make_path($dir,{mode=>0755}) or die "Can't create $dir: $!\n";
	}
}

open ERROR, ">", "$outdir/error.log";

###################################################################################################
## Run list_maker.pl
###################################################################################################

print ERROR "\n### list_maker.pl ###\n";

system("
	$path/list_maker.pl \\
	  --input @annot_files \\
	  --outdir $outdir \\
	  2>> $outdir/error.log
");

###################################################################################################
## Run get_homology.pl
###################################################################################################

my %linked_files;

### Checking the $outdir/PROT_SEQ for files
### unless the -prot command line switch is invoked

unless (@prot_files){

	opendir (FAA, "$outdir/PROT_SEQ") or die "\n\n[ERROR]\tCan't open $outdir/PROT_SEQ: $!\n\n";

	while (my $file = readdir(FAA)){
		if ($file =~ /\.faa$/){
			push (@prot_files, "$outdir/PROT_SEQ/$file");
		}
	}

}

link_files();
print ERROR "\n### get_homology.pl ###\n";

system("
	$path/get_homology.pl \\
	--input @prot_files \\
	--evalue $evalue \\
	--outdir $diamond_dir \\
	2>> $outdir/error.log
");

###################################################################################################
## Run get_synteny.pl
###################################################################################################

print ERROR "\n### get_synteny.pl ###\n";

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
					2>> $outdir/error.log
				");
			}
		}
		else{
			print "\t. . .\n";
		}
	}
}

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
");

###################################################################################################
## Create links files, karyotypes and nucleotide biases for Circos
###################################################################################################

## Karyotypes and nucleotide biases

print ERROR "\n### nucleotide_biases.pl ###\n";

my $ref = '';
if ($reference){
	$ref = "--reference $reference";
}

my $circos_plot = '';
if ($circos){
	$circos_plot = '--plot';
}

my $gap = $gaps[0];

system("
	$circos_path/nucleotide_biases.pl \\
	--outdir $outdir/CIRCOS \\
	--fasta $outdir/GENOME/*.fasta \\
	--winsize 10000 \\
	--step 5000 \\
	--circos \\
	--gap $gap \\
	$ref \\
	$circos_plot \\
	2>> $outdir/error.log
");

## links files

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
	");
}

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