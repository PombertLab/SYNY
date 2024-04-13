#!/usr/bin/perl
## Pombert Lab 2022

my $name = 'get_synteny.pl';
my $version = '0.7.3a';
my $updated = '2024-04-13';

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
		  -o RCCvsCCMP \\
		  -list_dir LISTS

OPTIONS:
*query and subject files must correspond to qseqid and sseqid respectively*
-ql (--query_list)		List (.list) file from query organism generated with list_maker.pl
-qb (--query_blast)		BLAST/DIAMOND homology searches in output format 6 for query_vs_subject
-sl (--subject_list)		List (.list) file from subject organism generated with list_maker.pl
-sb (--subject_blast)		BLAST/DIAMOND homology searches in output format 6 for subject_vs_query
-g (--gap)		Space allowed between adjacent genes [Default: 0]
-o (--outdir)		Output directory ## Writes detected gene pairs
-sd (--sumdir)		Summary output directory
--list_dir		Directory cotnaining .list files
EXIT

my $query_list;
my $query_blast;
my $subject_list;
my $subject_blast;
my $gap = 0;
my $outdir = 'SYNTENY';
my $sumdir = $outdir;
my $list_dir;

GetOptions(
	'ql|query_list=s' => \$query_list,
	'qb|query_blast=s' => \$query_blast,
	'sl|subject_list=s' => \$subject_list,
	'sb|subject_blast=s' => \$subject_blast,
	'g|gap=s' => \$gap,
	'o|outdir=s' => \$outdir,
	'sd|sumdir=s' => \$sumdir,
	'list_dir=s' => \$list_dir
);

my $pair_dir = $outdir."/PAIRS";
my $cluster_dir = $outdir."/CLUSTERS";

my @dirs = ($cluster_dir,$pair_dir,$sumdir);

foreach my $dir (@dirs){
	unless (-d $dir){
		make_path($dir,{mode=>0755});
	}
}

my $query_name = (fileparse($query_list,".list"))[0];
my $sub_name = (fileparse($subject_list,".list"))[0];

my $pair_file = ${query_name}."_vs_".${sub_name}.'.gap_'.$gap.'.pairs';
my $cluster_file = ${query_name}."_vs_".${sub_name}.'.gap_'.$gap.'.clusters';

###################################################################################################
## Parse DIAMOND files for homologs
###################################################################################################

### Loading locus_tag information from the query from its list file
my @loci = ();  ## Must keep track of the input order
                ## Sorting keys alphanumerically can break
                ## depending on the structure of the locus tags!

my $list_file = $list_dir.'/'.$query_name.'.list';
open LIST, '<', $list_file or die "Unable to read from $list_file: $!\n";
while (my $line = <LIST>){
	my @data = split("\t", $line);
	my $locus = $data[0];
	push(@loci, $locus);
}
close LIST;

### Loading diamond blast result for query
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

### Loading diamond blast result for subject
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

foreach my $q_locus (@loci){

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

			## Genes must be located on the same chromosome to be considered a neighbor
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

		## Different chromosomes; must be || not &&
		## => both queries and subject must be located on the same contigs as the previous ones
		if ($c_sc ne $p_sc || $c_qc ne $p_qc){
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
my @cluster_sizes;

foreach my $cluster (sort{$a <=> $b}(keys(%clusters))){

	my $cluster_num = sprintf("%0${padding}d",$cluster);
	my $ql_1 = (split("\t",$clusters{$cluster}[0]))[0];
	my $ql_2 = (split("\t",$clusters{$cluster}[-1]))[0];
	my $sl_1 = (split("\t",$clusters{$cluster}[0]))[2];
	my $sl_2 = (split("\t",$clusters{$cluster}[-1]))[2];

	print OUT "### Cluster $cluster_num; $ql_1 to $ql_2; $sl_1 to $sl_2; size = ".scalar(@{$clusters{$cluster}})." ###\n";
	push (@cluster_sizes, scalar(@{$clusters{$cluster}}));

	foreach my $line (@{$clusters{$cluster}}){
		print OUT $line."\n";
	}
}

close OUT;

###################################################################################################
## Generate summary file
###################################################################################################

my $summary_file = 'clusters_summary.tsv';
my $sum_outfile = $sumdir.'/'.$summary_file;
open SUM, '>>', $sum_outfile or die "Unable to write to $sum_outfile: $!\n";

my $total_cluster = scalar(@cluster_sizes);
my @len = sort({$b <=> $a}@cluster_sizes);

# median
my $median;
my $median_pos = $total_cluster/2;
if ($median_pos =~ /^\d+$/){
	$median = $len[$median_pos];
}
else {
	my $med1 = int($median_pos);
	my $med2 = $med1 + 1;
	$median = (($len[$med1] + $len[$med2])/2);
}
$median = sprintf("%.0f", $median);

# sum
my $sum;
foreach (@len){ $sum += $_; }

# longest/shortest
my $large = sprintf("%.0f", $len[0]);
my $small = sprintf("%.0f", $len[$#len]);

# average
my $average = sprintf("%.0f", ($sum/$total_cluster));

### N50, N75, N90
# thresholds to reach for N metrics
my $n50_td = $sum*0.5;
my $n75_td = $sum*0.75;
my $n90_td = $sum*0.9;
# sums to calculate
my $nsum50 = 0;
my $nsum75 = 0;
my $nsum90 = 0;
# metrics to capture
my $n50;
my $n75,
my $n90;

foreach (@len){
	$nsum50 += $_;
	if ($nsum50 >= $n50_td){
		$n50 = $_;
		last;
	}
}
foreach (@len){
	$nsum75 += $_;
	if ($nsum75 >= $n75_td){
		$n75 = $_;
		last
	}
}
foreach (@len){
	$nsum90 += $_;
	if ($nsum90 >= $n90_td){
		$n90 = $_;
		last;
	}
}

$n50 = sprintf ("%.0f", $n50);
$n75 = sprintf ("%.0f", $n75);
$n90 = sprintf ("%.0f", $n90);

print SUM '##### '.${query_name}."_vs_".${sub_name}.'; Gap = '.$gap.' #####'."\n";
print SUM '  Total number of proteins in clusters:'."\t".$sum."\n";
print SUM '  # of clusters:'."\t".$total_cluster."\n";
print SUM '  Longest:'."\t".$large."\n";
print SUM '  Shortest:'."\t".$small."\n";
print SUM '  Average cluster size:'."\t".$average."\n";
print SUM '  Median cluster size:'."\t".$median."\n";
print SUM '  N50:'."\t".$n50."\n";
print SUM '  N75:'."\t".$n75."\n";
print SUM '  N90:'."\t".$n90."\n";
print SUM "\n";

