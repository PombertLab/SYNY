#!/usr/bin/perl
# Pombert lab, 2022

my $name = 'run_syny.pl';
my $version = '0.5';
my $updated = '2022-06-18';

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

USAGE		${name} \\
		-a *.gbff \\
		-p *.prot \\
		-g 5 \\
		-o Chloropicon_synteny

OPTION
-a (--annot)	Annotation files (Supported files: gbff, gff, embl)
-p (--prot)	Protein files
-e (--evalue)	BLAST evalue cutoff [Default = 1e-10]
-f (--format)	Figure format(s) (.png,.eps,.jpg,.jpeg,.pdf,.pgf,.ps,.raw,.rgba,.svg,.svgz,.tif,.tiff) [Default: .svg]
-n (--name)	Figure file prefix [Default: PHYLOGENTIC_FIGURE]
-g (--gaps)	Allowable number of gaps between pairs [Default = 0]
-o (--outdir)	Output directory [Default = SYNY]
EXIT

die ("\n$usage\n") unless (@ARGV);

my @annot_files;
my @prot_files;
my $evalue = '1e-10';
my @formats;
my @gaps;
my $outdir = 'SYNY';

GetOptions(
	'a|annot=s@{1,}' => \@annot_files,
	'p|proteins=s@{1,}' => \@prot_files,
	'e|evalue=s' => \$evalue,
	'f|format=s{1,}' => \@formats,
	'g|gaps=s{0,}' => \@gaps,
	'o|outdir=s' => \$outdir,
);

unless(@gaps){
	@gaps = (0);
}

###################################################################################################
## Output directory creation and setup
###################################################################################################

my ($script,$path) = fileparse($0);

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
my $distance_dir = "$outdir/DIST_MAT";
my $figure_dir = "$outdir/FIGURES";

my @outdirs = ($list_dir,$prot_dir,$annot_dir,$diamond_dir,$db_dir,$synteny_dir,$distance_dir,$figure_dir);

foreach my $dir (@outdirs){
	unless (-d $dir){
		make_path($dir,{mode=>0755}) or die "Can't create $dir: $!\n";
	}
}

###################################################################################################
## Run list_maker.pl
###################################################################################################

system("
	$path/list_maker.pl \\
	  --input @annot_files \\
	  --out $list_dir
");

system("mv $list_dir/PROT_SEQ/*.faa $prot_dir; rm -r $list_dir/PROT_SEQ");
system("mv $list_dir/LISTS/*.list $list_dir; rm -r $list_dir/LISTS");
system("mv $list_dir/ANNOTATIONS/*.annotations $annot_dir; rm -r $list_dir/ANNOTATIONS");

###################################################################################################
## Run get_homology.pl
###################################################################################################

my %linked_files;
link_files();

system("
	$path/get_homology.pl \\
	--input $prot_dir/*.faa \\
	--evalue $evalue \\
	--outdir $diamond_dir
");

###################################################################################################
## Run get_synteny.pl
###################################################################################################

foreach my $annot_file_1 (sort(@annot_files)){
	my ($file_name_1,$dir,$ext) = fileparse($annot_file_1,'\..*');
	my $linked_file_1 = $linked_files{$annot_file_1};
	foreach my $annot_file_2 (sort(@annot_files)){
		if($annot_file_1 ne $annot_file_2){
			foreach my $gap (@gaps){
				my ($file_name_2,$dir,$ext) = fileparse($annot_file_2,'\..*');
				my $linked_file_2 = $linked_files{$annot_file_2};
				print "Identifying synteny between $file_name_1 and $file_name_2!\n";
				system("
					$path/get_synteny.pl \\
					--query_list $list_dir/$file_name_1.list \\
					--query_blast $diamond_dir/${file_name_1}_vs_${file_name_2}.diamond.6 \\
					--subject_list $list_dir/$file_name_2.list \\
					--subject_blast $diamond_dir/${file_name_2}_vs_${file_name_1}.diamond.6 \\
					--gap $gap \\
					--outdir $synteny_dir/gap_$gap
				");
			}
		}
	}
}

###################################################################################################
## Run id_conserved_regions.pl
###################################################################################################

system("
	$path/id_conserved_regions.pl \\
	--lists $list_dir \\
	--blasts $diamond_dir \\
	--outdir $conserved_dir
");

###################################################################################################
## Calculate distance matrices
###################################################################################################

foreach my $gap (@gaps){
	system("
		$path/make_distance_matrix.py \\
		-l $list_dir \\
		-n $synteny_dir/gap_$gap/NEIGHBORS \\
		-p $synteny_dir/gap_$gap/PAIRS \\
		-c $synteny_dir/gap_$gap/CLUSTERS \\
		-o $distance_dir/gap_$gap
	");
}

###################################################################################################
## Generate phylogenetic figures
###################################################################################################

foreach my $gap (@gaps){
	
	my $format = "";
	
	foreach my $form (@formats){
		$format .= "$form "
	}

	opendir(DIR,"$distance_dir/gap_$gap");
	
	foreach my $mat (readdir(DIR)){
		if (-f "$distance_dir/gap_$gap/$mat"){
			my ($prefix) = $mat =~ /(\w+)\.tsv/;
			system("
				$path/phylogenetic_figure_maker.py \\
				-d $distance_dir/gap_$gap/$mat \\
				-p $prefix \\
				-f $format \\
				-o $figure_dir/gap_$gap
			");
		}
	}

}

###################################################################################################
## Subroutines
###################################################################################################

sub link_files {

	opendir(DIR,$prot_dir);
	foreach my $file (readdir(DIR)){
		if(-f "$prot_dir/$file"){
			print "Pushing $file to prot_files\n";
			push(@prot_files,"$prot_dir/$file");
		}
	}

	if(scalar(@annot_files) != scalar(@prot_files)){
		die("The number of annotation files (".scalar(@annot_files).") does not equal the number of protein files (".scalar(@prot_files)."\n");
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