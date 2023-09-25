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
	  --out $list_dir \\
	  2>> $outdir/error.log
");

system("mv $list_dir/PROT_SEQ/*.faa $prot_dir; rm -r $list_dir/PROT_SEQ");
system("mv $list_dir/LISTS/*.list $list_dir; rm -r $list_dir/LISTS");
system("mv $list_dir/ANNOTATIONS/*.annotations $annot_dir; rm -r $list_dir/ANNOTATIONS");

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
			print $file."\n";
			push (@prot_files, $file)
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