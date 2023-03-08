#!/usr/bin/perl
#Pombert Lab, 2022

use strict; use warnings; use Getopt::Long qw(GetOptions); use File::Path qw(make_path); use File::Basename;

my $name = 'create_circos_links.pl';
my $version = '0.2a';
my $updated = '2022-02-07';

my $usage=<<"EXIT";
NAME		$name
VERSION		$version
UPDATED		$updated
SYNOPSIS	Creates links file for circos based on shared protein clusters

COMMAND		$name \\
			  RCC138.list \\
			  CCMP1205.list \\
			  RCCvsCCMP.gap0.clusters \\
			  links.txt
OPTIONS
-1 (--list_1)	1st list file for comparison (created by list_maker.pl)
-2 (--list_2)	2nd list file for comparison (created by list_maker.pl)
-c (--cluster)	Cluster file (created by get_synteny.pl)
-o (--outdir)	Output directory [Default = Syn_Circos]
EXIT

die "\t$usage\n" unless @ARGV;

print("Running $name\n");

my $list_1;
my $list_2;
my $cluster_file;
my $outdir = "Syn_Circos";

GetOptions(
	'1|list_1=s' => \$list_1,
	'2|list_2=s' => \$list_2,
	'c|cluster=s' => \$cluster_file,
	'o|outdir=s' => \$outdir,
);

unless(-d $outdir){
	make_path($outdir,{mode=>0755}) or die "Unable to make directory $outdir: $!\n";
}

open RCC, "<", "$list_1" or die("Unable to open list file $list_1: $!\n");
my %RCC; my %CCMP;
while (my $line = <RCC>){
	chomp $line;
	my @cols = split("\t", $line);
	foreach (@cols){
		push ( @{$RCC{$cols[0]}}, $_);
		# [0] => locus_tag, [1] => chromosome
		# [2] => start, [3] => end
		# [4] => strand, [5] => gene number
	}
}

open CCMP, "<", "$list_2" or die("Unable to open list file $list_2: $!\n");
while (my $line = <CCMP>){
	chomp $line;
	my @cols = split("\t", $line);
	$cols[1] =~ s/\.1//;
	foreach (@cols){
		push ( @{$CCMP{$cols[0]}}, $_);
	}
}

my ($file_name_1,$dir_1,$ext_1) = fileparse($list_1,'\..*');
my ($file_name_2,$dir_2,$ext_2) = fileparse($list_2,'\..*');
my $outfile = "${file_name_1}_vs_${file_name_2}.txt";
open OUT, ">", "$outdir/$outfile";
print OUT '#species1 start end species2 start end'."\n";

open CLU, "<", "$cluster_file";
while (my $line = <CLU>){
	chomp $line;
	if ($line =~ /### Cluster (\d+); (\w+) to (\w+); (\S+) to (\S+);/){
		my $cluster = $1;
		my $rloc1 = $2;
		my $rloc2 = $3;
		my $cloc1 = $4;
		my $cloc2 = $5;

		print OUT "$RCC{$rloc1}[1] $RCC{$rloc2}[2] $RCC{$rloc1}[3] ";
		print OUT "$CCMP{$cloc1}[1] $CCMP{$cloc2}[2] $CCMP{$cloc1}[3] ";
		print OUT '### Cluster '."$cluster\n";
	}
}
