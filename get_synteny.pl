#!/usr/bin/perl
## Pombert Lab 2022

my $name = 'get_synteny.pl';
my $version = '0.7.2';
my $updated = '2022-06-22';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Path qw(make_path);
use File::Basename;

my $usage = <<"EXIT";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Detects pairs of genes that are shared between genomes based on BLAST/DIAMOND homology searches,
		then reconstructs clusters based on these pairs.

USAGE		$name \\
		  -ql RCC138.list \\
		  -qb RCC_vs_CCMP.diamond.blastp.6 \\
		  -sl CCMP1205.list \\
		  -sb CCMP_vs_RCC.diamond.blastp.6
		  -gap 10 \\
		  -o RCCvsCCMP

OPTIONS:
*query and subject files must correspond to qseqid and sseqid respectively*
-ql (--query_list)		List (.list) file from query organism generated with list_maker.pl
-qb (--query_blast)		BLAST/DIAMOND homology searches in output format 6 for query_vs_subject
-sl (--subject_list)		List (.list) file from subject organism generated with list_maker.pl
-sb (--subject_blast)		BLAST/DIAMOND homology searches in output format 6 for subject_vs_query
-g (--gap)		Space allowed between adjacent genes [Default: 0]
-o (--outdir)		Output directory ## Writes detected gene pairs
EXIT

my $query_list;
my $query_blast;
my $subject_list;
my $subject_blast;
my $gap = 0;
my $outdir = 'SYNTENY';

GetOptions(
	'ql|query_list=s' => \$query_list,
	'qb|query_blast=s' => \$query_blast,
	'sl|subject_list=s' => \$subject_list,
	'sb|subject_blast=s' => \$subject_blast,
	'g|gap=s' => \$gap,
	'o|outdir=s' => \$outdir,
);

my $pair_dir = $outdir."/PAIRS";
my $cluster_dir = $outdir."/CLUSTERS";

my @dirs = ($cluster_dir,$pair_dir);

foreach my $dir (@dirs){
	unless (-d $dir){
		make_path($dir,{mode=>0755});
	}
}

my $query_name = (fileparse($query_list,".list"))[0];
my $sub_name = (fileparse($subject_list,".list"))[0];

my $pair_file = ${query_name}."_vs_".${sub_name}.".pairs";
my $cluster_file = ${query_name}."_vs_".${sub_name}.".clusters";

###################################################################################################
## Parse DIAMOND files for homologs
###################################################################################################

open QB, "<", $query_blast or die "Unable to read from $query_blast: $!\n";
my %q_blast_hits;
while (my $line = <QB>){
	chomp($line);
	my @data = split("\t",$line);
	my $q_locus = $data[0];
	my $s_locus = $data[1];
	my $e_value = $data[10];
	unless($q_blast_hits{$s_locus}){
		@{$q_blast_hits{$s_locus}} = ($q_locus,$e_value);
	}
	elsif($e_value < @{$q_blast_hits{$s_locus}}[1]){
		@{$q_blast_hits{$s_locus}} = ($q_locus,$e_value);
	}
}
close QB;

open SB, "<", $subject_blast or die "Unable to read from $subject_blast: $!\n";
my %s_blast_hits;
while (my $line = <SB>){
	chomp($line);
	my @data = split("\t",$line);
	my $s_locus = $data[0];
	my $q_locus = $data[1];
	my $e_value = $data[10];
	unless($s_blast_hits{$q_locus}){
		@{$s_blast_hits{$q_locus}} = ($s_locus,$e_value);
	}
	elsif($e_value < @{$s_blast_hits{$q_locus}}[1]){
		@{$s_blast_hits{$q_locus}} = ($s_locus,$e_value);
	}
}
close QB;

###################################################################################################
## Verfiy bidirectionality of DIAMOND results
###################################################################################################

my %q_bd_hits;
foreach my $q_s_locus (keys(%q_blast_hits)){
	my $q_q_locus = @{$q_blast_hits{$q_s_locus}}[0];
	if($s_blast_hits{$q_q_locus}){
		my $s_s_locus = @{$s_blast_hits{$q_q_locus}}[0];
		if ($s_s_locus eq $q_s_locus){
			$q_bd_hits{$q_q_locus} = $s_s_locus;
		}
	}
}

###################################################################################################
## Store query and subject info locally
###################################################################################################

open QL, "<", $query_list or die "Unable to read from $query_list: $!\n";

my %query_info;
my $counter = 1;
while (my $line = <QL>){
	chomp($line);
	my ($locus,$accession,$start,$end,$strand,$number) = split("\t",$line);
	@{$query_info{$locus}} = ($accession,$strand,$counter);
	$counter += 1;
}

close QL;

open SL, "<", $subject_list or die "Unable to read from $subject_list: $!\n";

my %subject_info;
$counter += 1;
while (my $line = <SL>){
	chomp($line);
	my ($locus,$accession,$start,$end,$strand,$number) = split("\t",$line);
	@{$subject_info{$locus}} = ($accession,$strand,$counter);
	$counter += 1;
}

close SL;

###################################################################################################
## Generate pairs file
###################################################################################################

## Previous query locus
my $p_ql;
## Previous query chromosome
my $p_qc;
## Previous query strand
my $p_qs;
## Previous query protein number
my $p_qn;

## Previous subject locus
my $p_sl;
## Previous subject chromosome
my $p_sc;
## Previous subject strand
my $p_ss;
## Previous subject protein number
my $p_sn;

open OUT, ">", $pair_dir."/".$pair_file or die "Unable to read from $pair_dir/$pair_file: $!\n";

foreach my $q_locus (sort(keys(%q_bd_hits))){
	
	## Current query locus
	my $c_ql = $q_locus;
	## Current query chromosome
	my $c_qc = $query_info{$c_ql}[0];
	## Current query strand
	my $c_qs = $query_info{$c_ql}[1];
	## Current query protein number
	my $c_qn = $query_info{$c_ql}[2];
	
	## Current subject locus
	my $c_sl;
	## Current subject chromosome
	my $c_sc;
	## Current subject strand
	my $c_ss;
	## Current subject protein number
	my $c_sn;

	## Gene must have a BLAST hit to be considered a neighbor
	if ($q_bd_hits{$c_ql}){

		$c_sl = $q_bd_hits{$c_ql};
		$c_sc = $subject_info{$c_sl}[0];
		$c_ss = $subject_info{$c_sl}[1];
		$c_sn = $subject_info{$c_sl}[2];

		## Both a query and subject locus is required
		if($p_ql && $p_sl){

			## Genes must be located on the same chromomsome to be considered a nieghbor
			if ($c_qc eq $p_qc){
				if ($c_sc eq $p_sc){
					

					if ((((($p_qn - $c_qn)**2)**(1/2))-1) <= $gap ){
						if ((((($p_sn - $c_sn)**2)**(1/2))-1) <= $gap ){
							print OUT $p_ql."\t".$p_qn."\t";
							print OUT $c_ql."\t".$c_qn."\t";
							print OUT $c_qc."\t";
							print OUT $p_sl."\t".$p_sn."\t";
							print OUT $c_sl."\t".$c_sn."\t";
							print OUT $c_sc."\t";
							my $q_o = $p_qs.$c_qs;
							my $s_o = $p_ss.$c_ss;
							print OUT $q_o."/".$s_o."\n";
						}
					}
				}
			}
		}
	}

	$p_ql = $c_ql;
	$p_qc = $c_qc;
	$p_qs = $c_qs;
	$p_qn = $c_qn;
	
	$p_sl = $c_sl;
	$p_sc = $c_sc;
	$p_ss = $c_ss;
	$p_sn = $c_sn;
}

close OUT;


###################################################################################################
## Generate clusters file
###################################################################################################

open IN, "<", $pair_dir."/".$pair_file or die "Unable to read from $pair_dir/$pair_file: $!\n";

undef $p_ql;
undef $p_qn;
undef $p_qs;
undef $p_qc;

undef $p_sn;
undef $p_sn;
undef $p_ss;
undef $p_sc;

my %clusters;
my $cluster_number = 1;
while (my $line = <IN>){
	chomp($line);
	my @data = split("\t",$line);
	
	## Current query 1 locus
	my $c_ql_1 = $data[0];
	## Current query 1 number
	my $c_ql_1_n = $data[1];
	## Current query strand
	my $c_ql_1_s = substr($data[10],0,1);
	## Current query 2 locus
	my $c_ql_2 = $data[2];
	## Current query 2 number
	my $c_ql_2_n = $data[3];
	## Current query 2
	my $c_ql_2_s = substr($data[10],1,1);
	## Current query chromosome
	my $c_qc = $data[4];

	my $c_sl_1 = $data[5];
	my $c_sl_1_n = $data[6];
	my $c_sl_1_s = substr($data[10],3,1);
	my $c_sl_2 = $data[7];
	my $c_sl_2_n = $data[8];
	my $c_sl_2_s = substr($data[10],4,1);
	my $c_sc = $data[9];

	## First locus will not have access to previous locus information
	if ($p_ql){

		## The gap between the previous query and current query is abs() -1
		my $query_gap = ((($c_ql_1_n - $p_qn)**2)**(1/2)) - 1;
		## The gap between the previous subject and current subject is abs() -1
		my $sub_gap = ((($c_sl_1_n - $p_sn)**2)**(1/2)) - 1;

		## Different chromosomes
		if ($c_sc ne $p_sc && $c_qc ne $p_qc){
			$cluster_number ++;
			push (@{$clusters{$cluster_number}},$c_ql_1."\t".$c_ql_1_s."\t".$c_sl_1."\t".$c_sl_1_s);
			push (@{$clusters{$cluster_number}},$c_ql_2."\t".$c_ql_2_s."\t".$c_sl_2."\t".$c_sl_2_s);
		}
		## No gaps
		elsif ($query_gap == -1 && $sub_gap == -1){
			push (@{$clusters{$cluster_number}},$c_ql_2."\t".$c_ql_2_s."\t".$c_sl_2."\t".$c_sl_2_s);
		}
		## Acceptable query gap size
		elsif ( $query_gap <= $gap){
			## Acceptable subject gap size
			if ( $sub_gap <= $gap){
				push (@{$clusters{$cluster_number}},"## Gap of ".$query_gap."\t\t"."Gap of ".$sub_gap." ##\t");
				push (@{$clusters{$cluster_number}},$c_ql_1."\t".$c_ql_1_s."\t".$c_sl_1."\t".$c_sl_1_s);
				push (@{$clusters{$cluster_number}},$c_ql_2."\t".$c_ql_2_s."\t".$c_sl_2."\t".$c_sl_2_s);
			}
			## Unacceptable subject gap size
			else {
				$cluster_number ++;
				push (@{$clusters{$cluster_number}},$c_ql_1."\t".$c_ql_1_s."\t".$c_sl_1."\t".$c_sl_1_s);
				push (@{$clusters{$cluster_number}},$c_ql_2."\t".$c_ql_2_s."\t".$c_sl_2."\t".$c_sl_2_s);
			}
		}
		## Unacceptable query gap size
		else {
			$cluster_number ++;
			push (@{$clusters{$cluster_number}},$c_ql_1."\t".$c_ql_1_s."\t".$c_sl_1."\t".$c_sl_1_s);
			push (@{$clusters{$cluster_number}},$c_ql_2."\t".$c_ql_2_s."\t".$c_sl_2."\t".$c_sl_2_s);
		}
	}
	else {
		push (@{$clusters{$cluster_number}},$c_ql_1."\t".$c_ql_1_s."\t".$c_sl_1."\t".$c_sl_1_s);
		push (@{$clusters{$cluster_number}},$c_ql_2."\t".$c_ql_2_s."\t".$c_sl_2."\t".$c_sl_2_s);
	}

	$p_ql = $c_ql_2;
	$p_qn = $c_ql_2_n;
	$p_qs = $c_ql_2_s;
	$p_qc = $c_qc;

	$p_sl = $c_sl_2;
	$p_sn = $c_sl_2_n;
	$p_ss = $c_sl_2_s;
	$p_sc = $c_sc;
}

close IN;

open OUT, ">", $cluster_dir."/".$cluster_file or die "Unable to write to $cluster_dir/$cluster_file: $!\n";

my $padding = $cluster_number =~ s/\d//g;
$padding ++;
foreach my $cluster (sort{$a <=> $b}(keys(%clusters))){
	my $cluster_num = sprintf("%0${padding}d",$cluster);
	my $ql_1 = (split("\t",$clusters{$cluster}[0]))[0];
	my $ql_2 = (split("\t",$clusters{$cluster}[-1]))[0];
	my $sl_1 = (split("\t",$clusters{$cluster}[0]))[2];
	my $sl_2 = (split("\t",$clusters{$cluster}[-1]))[2];
	print OUT "### Cluster $cluster_num; $ql_1 to $ql_2; $sl_1 to $sl_2; size = ".scalar(@{$clusters{$cluster}})." ###\n";
	foreach my $line (@{$clusters{$cluster}}){
		print OUT $line."\n";
	}
}

close OUT;