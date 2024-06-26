#!/usr/bin/env perl
## Pombert Lab, 2022

my $name = "get_homology.pl";
my $version = "0.1.8a";
my $updated = "2024-05-27";

use warnings;
use strict;
use Getopt::Long qw(GetOptions);
use File::Basename;
use File::Path qw(make_path);

#########################################################################
### Command line options
#########################################################################

my $usage = << "EXIT";
NAME        ${name}
VERSION     ${version}
UPDATED     ${updated}
SYNOPSIS    Creates DIAMOND databases based on input protein files and performs round-robin BLASTP homology searches

USAGE       ${name} \\
              -i *.prot \\
              -a *.annotations \\
              -e 1e-40 \\
              -o DIAMOND

OPTIONS
-i (--input)    Input protein files
-l (--list)     Directory containing tab-delimited lists of features [Default = LISTS]
-e (--evalue)   Evalue [Default = 1e-10]
-o (--outdir)   Diamond searches output directory [Default = DIAMOND]
-s (--shared)   Shared proteins output directory [Default = SHARED]
-t (--threads)  Number of threads to use [Default = 8]
-v (--version)  Show script version
EXIT

unless (@ARGV){
	print "\n$usage\n";
	exit(0);
};

my @input_files;
my $e_value = "1e-10";
my $outdir = "DIAMOND";
my $shared_dir = "SHARED";
my $list_dir = "LISTS";
my $threads = 8;
my $sc_version;
GetOptions(
	'i|input=s@{2,}' => \@input_files,
	'e|evalue=s' => \$e_value,
	'o|outdir=s' => \$outdir,
	's|shared=s' => \$shared_dir,
	'l|list=s' => \$list_dir,
	't|threads=i' => \$threads,
	'v|version' => \$sc_version
);

#########################################################################
### Version
#########################################################################

if ($sc_version){
    print "\n";
    print "Script:     $name\n";
    print "Version:    $version\n";
    print "Updated:    $updated\n\n";
    exit(0);
}

#########################################################################
### Req. check
#########################################################################

my $diamond_check = `echo \$(command -v diamond)`;
chomp $diamond_check;
if ($diamond_check eq ''){
	print STDERR "\n[E]: Cannot find diamond. Please install diamond in your \$PATH. Exiting..\n\n";
	exit;
}

#########################################################################
### Output dir/subdirs
#########################################################################

unless (-d $outdir){
	make_path($outdir,{mode=>0755});
}

my $db_dir = "$outdir/DB";
unless (-d $db_dir){
	make_path($db_dir,{mode=>0755});
}

unless (-d $shared_dir){
	make_path($shared_dir,{mode=>0755});
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
			  --threads $threads \\
			  --in $file \\
			  --db $db_dir/$file_prefix \\
			  --quiet
		") == 0 or checksig();
	}
	$db_files{$file_prefix} = "$db_dir/$file_prefix";
}

###################################################################################################
# Running homology searches
###################################################################################################

my @prefixes;
my @diamond_results;

foreach my $file (sort(@input_files)){

	my ($file_name, $path, $suffix) = fileparse($file);
	my $file_prefix = (split('\.',$file_name))[0];
	push (@prefixes, $file_prefix);

	print "\nRunning DIAMOND BLASTP on $file\n";

	foreach my $db_prefix (keys(%db_files)){

		if($db_prefix ne $file_prefix){

			my $blast_file = "$outdir/${file_prefix}_vs_${db_prefix}.diamond.6";
			push (@diamond_results, $blast_file);
			
			print "\tusing DB $db_prefix\n";

			unless(-f $blast_file){
				system ("diamond blastp \\
						--threads $threads \\
						-d $db_files{$db_prefix} \\
						-q $file \\
						-o $blast_file \\
						-e $e_value \\
						-k 1 \\
						-f 6 \\
						--quiet
				") == 0 or checksig();
			}
			
		}
		else{
			print "\t. . .\n";
		}
	}
}

###################################################################################################
# Parsing homology searches / getting shared and unique proteins based on e-value cutoff
###################################################################################################

##### Parsing diamond output files

my %hom_results;

foreach my $diamond (@diamond_results){

	open DD, '<', $diamond or die "Can't open $diamond: $!\n";
	my ($basename) = fileparse($diamond);
	$basename =~ s/\.diamond\.6$//;

	while (my $line = <DD>){

		chomp $line;

		if ($line =~ /^#/){
			next;
		}
		else {

			my @data = split ("\t", $line);

			my $query = $data[0];
			my $subject = $data[1];
			my $eval = $data[10];

			if ($eval <= $e_value){
				
				if (exists $hom_results{$basename}{$query}){
					## Keeping only the best evalues
					if ($eval < $hom_results{$basename}{$query}{'evalue'}){
						$hom_results{$basename}{$query}{'subject'} = $subject;
						$hom_results{$basename}{$query}{'evalue'} = $eval;
					}
				}

				else {
					$hom_results{$basename}{$query}{'subject'} = $subject;
					$hom_results{$basename}{$query}{'evalue'} = $eval;
				}

			}

		}

	}

	close DD;

}

##### Creating lists for each species used

foreach my $prefix (@prefixes){

	my $protein_list = $list_dir.'/'.$prefix.'.list';

	open LIST, '<', $protein_list or die "Can't open $protein_list: $!\n";

	my $uniques = $shared_dir.'/'.$prefix.'.uniques.tsv';
	my $shared = $shared_dir.'/'.$prefix.'.shared.tsv';
	my $all = $shared_dir.'/'.$prefix.'.all.tsv';

	open UNI, '>', $uniques or die "Can't create $uniques: $!\n"; 
	open SHA, '>', $shared or die "Can't create $shared: $!\n";
	open ALL, '>', $all or die "Can't create $all: $!\n";

	## Printing header
	for my $fh (\*SHA,\*ALL){
		print $fh '# '.$prefix."\t".'Description';
		foreach my $subject (sort @prefixes){
			if ($subject ne $prefix){
				print $fh "\t".$subject.' locus';
				print $fh "\t".$subject.' evalue';
			}
		}
		print $fh "\n";
	}

	# Working on data
	while (my $line = <LIST>){

		chomp $line;

		my @data = split("\t", $line);
		my $protein = $data[0]; 
		my $prot_desc = $data[-1];

		my $found;
		my $line_out = $protein."\t".$prot_desc;

		foreach my $key (sort (keys %hom_results)){

			## Print only relevant entrys from the master results database
			if ($key =~ /^$prefix/){

				if (exists $hom_results{$key}{$protein}){
					$found = 1;
					$line_out .= "\t$hom_results{$key}{$protein}{subject}";
					$line_out .= "\t$hom_results{$key}{$protein}{evalue}";
				}
				else {
					$line_out .= "\t---" x 2;
				}

			}

		}

		if ($found){
			print SHA $line_out."\n";
			print ALL $line_out."\n";
		}
		else {
			print UNI $protein."\t".$prot_desc."\n";
			print ALL $line_out."\n";
		}

		$found = undef;

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