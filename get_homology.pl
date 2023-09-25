#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);
use File::Basename;
use File::Path qw(make_path);

my $name = "get_homology.pl";
my $version = "0.1.4";
my $updated = "2023-02-03";
my $usage = << "EXIT";
NAME		$name
VERSION		$version
UPDATED		$updated
SYNOPSIS	Creates DIAMOND databases based on input protein files and performs round-robin BLASTP homology searches

USAGE		$name \\
		  -i *.prot \\
		  -e 1e-40 \\
		  -o DIAMOND

OPTIONS
-i (--input)	Input protein files
-e (--evalue)	Evalue [Default = 1e-10]
-o (--outdir)	Output directory [Default = DIAMOND]
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
	make_path($db_dir,{mode=>0755});
}

###################################################################################################
# Creating diamond databases
###################################################################################################

print "\n";

my @db_files;
print "Creating DIAMOND DB\n";
foreach my $file (@input_files){

	my ($basename, $path, $suffix) = fileparse($file);
	my ($file_name) = $basename =~ /^(\S+)\.(fasta|faa|fa|prot)$/;

	unless (-f "$db_dir/$file_name.dmnd"){
		print "\t$file_name\n";
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

	my ($basename, $path, $suffix) = fileparse($file);
	my ($file_name) = $basename =~ /^(\S+)\.(fasta|faa|fa|prot)$/;

	print "\nRunning BLASTP on $file\n";

	foreach my $db_file (sort(@db_files)){

		my ($db_name, @aux) = fileparse($db_file,'\..*');

		if($db_name ne $file_name){

			my $blast_file = "$outdir/${file_name}_vs_${db_name}.diamond.6";
			
			print "\tusing DB $db_name\n";

			unless(-f $blast_file){
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
		else{
			print "\t. . .\n";
		}
	}
}