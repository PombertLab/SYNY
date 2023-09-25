#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);
use File::Basename;
use File::Path qw(make_path);

my $name = "get_homology.pl";
my $version = "0.1.5";
my $updated = "2023-09-25";
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

my %db_files;
print "Creating DIAMOND DB\n";
foreach my $file (@input_files){

	my ($file_name, $path, $suffix) = fileparse($file);

	my $file_prefix = (split('\.',$file_name))[0];

	unless (-f "$db_dir/$file_prefix.dmnd"){
		print "\t$file_name\n";
		system ("
			diamond makedb \\
			  --in $file \\
			  --db $db_dir/$file_prefix &> /dev/null
		");
	}
	$db_files{$file_prefix} = "$db_dir/$file_prefix";
}

###################################################################################################
# Running homology searches
###################################################################################################

foreach my $file (sort(@input_files)){

	my ($file_name, $path, $suffix) = fileparse($file);
	my $file_prefix = (split('\.',$file_name))[0];

	print "\nRunning BLASTP on $file\n";

	foreach my $db_prefix (keys(%db_files)){

		if($db_prefix ne $file_prefix){

			my $blast_file = "$outdir/${file_prefix}_vs_${db_prefix}.diamond.6";
			
			print "\tusing DB $db_prefix\n";

			unless(-f $blast_file){
				system ("diamond blastp \\
						-d $db_files{$db_prefix} \\
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