#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);
use File::Basename;
use File::Path qw(make_path);

my $name = "get_homology.pl";
my $version = "0.1b";
my $updated = "2022-02-14";
my $usage = << "EXIT";
NAME		$name
VERSION		$version
UPDATED		$updated
SYNOPSIS	The purpose of this script is to create DIAMOND databases based on input protein files
		and then cross blastp all files

USAGE		$name \\
		  -i *.prot \\
		  -e 1e-40 \\
		  -o DIAMOND

OPTIONS
-i (--input)	Input protein files
-e (--evalue)	Evalue [Default = 1e-10]
-o (--outidr)	Output directory [Default = DIAMOND]
EXIT

die("\n$usage\n") unless (@ARGV);

###################################################################################################
# Setting up script variables
###################################################################################################

my @input_files;
my $e_value = "1e-10";
my $outdir = "DIAMOND";

GetOptions(
	'i|input=s@{2,}' => \@input_files,
	'e|evalue=s' => \$e_value,
	'o|outdir=s' => \$outdir,
);

my $diamond_version = `which diamond`;
if($diamond_version eq ''){
	die("Make sure diamond is installed and in the PATH\n");
}

unless(-d $outdir){
	make_path($outdir,{mode=>0755});
}

my $db_dir = "$outdir/DB";
unless(-d $db_dir){
	mkdir(0755,$db_dir);
}

###################################################################################################
# Creating diamond databases
###################################################################################################

my @db_files;
foreach my $file (@input_files){
	my ($file_name, @aux) = fileparse($file,'\..*');

	unless (-f "$db_dir/$file_name.dmnd"){
		print("Creating DIAMOND DB $db_dir/$file_name\n");
		system ("
			diamond makedb \\
			  --in $file \\
			  --db $db_dir/$file_name &> /dev/null
		");
	}
	push(@db_files,"$db_dir/$file_name");
}

###################################################################################################
# Running homology searches
###################################################################################################

foreach my $file (sort(@input_files)){

	my ($file_name, @aux) = fileparse($file,'\..*');

	foreach my $db_file (sort(@db_files)){

		my ($db_name, @aux) = fileparse($db_file,'\..*');

		if($db_name ne $file_name){

			my $blast_file = "$outdir/${file_name}_vs_${db_name}.diamond.6";
			
			unless(-f $blast_file){
				print("Running BLASTP on $file using $db_file DB\n");
				system ("diamond blastp \\
						-d $db_file \\
						-q $file \\
						-o $blast_file \\
						-e $e_value \\
						-k 1 \\
						-f 6 &> /dev/null
				")
			}
			
		}
	}
}