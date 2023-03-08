#!/usr/bin/perl

use warnings; use strict; use Getopt::Long qw(GetOptions); use File::Basename; use File::Path qw(make_path); use Cwd qw(abs_path);

my $name = "synteny.pl";
my $version = "0.0a";
my $updated = "2022/01/29";
my $usage = << "EXIT";
NO INPUT PROVIDED
EXIT

die("\n$usage\n") unless (@ARGV);

###################################################################################################
# Setting up script variables
###################################################################################################

my @input_files;
my $outdir = "SYNTENY";

GetOptions(
	'i|input=s@{2,}' => \@input_files,
	'o|outdir=s' => \$outdir,
);

die("synteny.pl requires at least 2 files for comparisons.\n") unless(scalar(@input_files) > 1);

my $diamond_version = `which diamond`;
if($diamond_version eq ''){
	die("Make sure diamond is installed and in the PATH\n");
}

unless(-d $outdir){
	make_path(0755,$outdir) or die("Could not create output directory: $!\n");
}

my $temp_dir = "$outdir/_temp_";
my $blast_dir = "$temp_dir/diamond";
my $db_dir = "$blast_dir/dbs";
unless(-d $db_dir){
	make_path(0755,$db_dir);
}

###################################################################################################
# Checking file extension for any improper inputs
###################################################################################################

my @prot_files;
foreach my $file (@input_files){
	my ($file_name,$dir,$ext) = fileparse($file,'\..*');
	if($ext eq ".prot"){
		push(@prot_files,abs_path($file));
		next;
	}
	else{
		## Run transformation file
	}
}

###################################################################################################
# Creating diamond databases
###################################################################################################

my @db_files;
my %locus_to_order;
my %order_to_locus;
foreach my $file (@prot_files){
	open IN, "<", $file;
	my ($file_name, @aux) = fileparse($file,'\..*');
	
	my @tags;
	while (my $line = <IN>){
		chomp($line);
		if ($line =~ />(\w+)/){
			push(@tags,$1);
		}
	}

	my $counter = 1;
	foreach my $tag (sort(@tags)){
		$locus_to_order{$file_name}{$tag} = $counter;
		$order_to_locus{$file_name}{$counter} = $tag;
		$counter++;
	}

	unless (-f "$db_dir/$file_name.dmnd"){
		system ("diamond makedb \\
				--in $file \\
				--db $db_dir/$file_name
		");
	}
	push(@db_files,"$db_dir/$file_name");
}

###################################################################################################
# Running homology searches
###################################################################################################

my %synteny;
foreach my $file (sort(@prot_files)){

	my ($file_name, @aux) = fileparse($file,'\..*');

	foreach my $db_file (sort(@db_files)){

		my ($db_name, @aux) = fileparse($db_file,'\..*');

		if($db_name ne $file_name){

			my $blast_file = "$blast_dir/${file_name}_vs_${db_name}.diamond.6";
			
			unless(-f $blast_file){
				system ("diamond blastp \\
						-d $db_file \\
						-q $file \\
						-o $blast_file \\
						-k 1 \\
						-f 6
				")
			}
			
			open BLAST ,"<", $blast_file;
			open SYNTENY, ">", "$temp_dir/${file_name}_vs_${db_name}.synteny";

			my $previous_query_locus;
			my $previous_db_locus;
			my $locus_tracker = 1;
			my $cluster_length = 1;
			my $cluster_start;
			my $cluster_end;
			
			while (my $line = <BLAST>){
				chomp($line);
				RETURN:
				my ($query_locus,$db_locus,@junk) = split("\t",$line);

				## Cluster identification
				if($previous_db_locus){

					my $prev_index = $locus_to_order{$db_name}{$previous_db_locus}+1;
					my $curr_index = $locus_to_order{$db_name}{$db_locus};
					if($locus_to_order{$db_name}{$db_locus} == $locus_to_order{$db_name}{$previous_db_locus}+1){
						$cluster_start = $previous_query_locus;
						$cluster_length++;
					}
					elsif($cluster_length > 1){
						$cluster_end = $previous_query_locus;
						# print("$cluster_length sequential genes in cluster from $cluster_start to $cluster_end\n");
						$cluster_length = 1;
					}
				}

				## Homology matches
				if($order_to_locus{$file_name}{$locus_tracker} eq $query_locus){
					$synteny{$file_name}{$query_locus} .= "\t$db_locus";
					print SYNTENY ("$query_locus\t$db_locus\n");
					$locus_tracker++;
				}
				else{
					$synteny{$file_name}{$order_to_locus{$file_name}{$locus_tracker}} .= tab($db_locus);
					print SYNTENY ($order_to_locus{$file_name}{$locus_tracker}."\n");
					$locus_tracker++;
					goto RETURN;
				}
				$previous_db_locus = $db_locus;
				$previous_query_locus = $query_locus;
			}

			close BLAST;
			close SYNTENY;

			open SYNTENY, ">", "$temp_dir/${file_name}_all.synteny";

			foreach my $locus (sort(keys(%{$synteny{$file_name}}))){
				my $synts = $synteny{$file_name}{$locus};
				print SYNTENY ("$locus$synts\n");
			}
		}
	}
}

print("\n");

sub tab {
	my $val = length($_[0]);
	my $tabs = int($val/4);

	if(($val/4) - int($val/4) > 0){
		$tabs++;
	}

	my $tab = "";
	for (my $i = 0; $i < $tabs; $i++){
		$tab .= "\t";
	}
	return ($tab);
}
