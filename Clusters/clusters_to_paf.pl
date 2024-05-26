#!/usr/bin/env perl
## Pombert Lab, 2024

my $name = 'clusters_to_paf.pl';
my $version = '0.2';
my $updated = '2024-04-22';

use strict;
use warnings;
use File::Basename;
use File::Path qw(make_path);
use Getopt::Long qw(GetOptions);

my $usage = << "USAGE";
NAME        ${name}
VERSION     ${version}
UDPATED     ${updated}
SYNOPSIS    Generates pseudo-PAF files from clusters found with SYNY

COMMAND     ${name} \\
              --lists *.list \\
              --clusters *.clusters \\
              --fasta *.fasta

OPTIONS:
-l (--lists)        Feature list files generated by SYNY
-c (--clusters)     Clusters files generated by SYNY
-f (--fasta)        FASTA files
-o (--outdir)       Output directory [Default: ./PAF]
USAGE

unless (@ARGV){
	print "\n$usage\n";
	exit(0);
};

my @clusters;
my @fastas;
my @lists;
my $outdir = './PAF';
GetOptions(
    'c|clusters=s@{1,}' => \@clusters,
    'f|fasta=s@{1,}' => \@fastas,
    'l|lists=s@{1,}' => \@lists,
    'o|outdir=s' => \$outdir
);

###################################################################################################
## Output dir/subdir
###################################################################################################

unless (-d $outdir){
    make_path($outdir,{mode=>0755}) or die("Can't create output directory $outdir: $!\n");
}

###################################################################################################
## Loading fasta and annotations
###################################################################################################

my %annotations;
my %lengths;

## Grabbing sequence lengths from FASTA files
while (my $fasta = shift@fastas){

    open FASTA, '<', $fasta or die "Can't open $fasta: $!\n";

    my %sequences;
    my $seqname;

    while (my $line = <FASTA>){
        chomp $line;
        if ($line =~ /^>(\S+)/){
            $seqname = $1;
        }
        else {
            $sequences{$seqname} .= $line;
        }
    }

    for my $seq (keys %sequences){
        my $len = length $sequences{$seq};
        $lengths{$seq} = $len;
    }

}

## Grabbing data from feature lists
while (my $list = shift@lists){

    open LIST, '<', $list or die "Can't open $list: $!\n";

    while (my $line = <LIST>){

        chomp $line;

        my @data = split("\t", $line);
        my $locus_tag = $data[0];
        my $contig = $data[1];
        my $start = $data[2];
        my $end = $data[3];
        my $strand = $data[4];
        my $gene_num = $data[5];
        my $description = $data[6];

        $annotations{$locus_tag}{'contig'} = $contig;
        $annotations{$locus_tag}{'start'} = $start;
        $annotations{$locus_tag}{'end'} = $end;
        $annotations{$locus_tag}{'strand'} = $strand;
        $annotations{$locus_tag}{'gene_num'} = $gene_num;
        $annotations{$locus_tag}{'description'} = $description;

    }

}


###################################################################################################
## Converting clusters to PAF
###################################################################################################

my $map_quality = 255; ## 255 = missing

while (my $cluster = shift @clusters){

    my ($basename,$path) = fileparse($cluster);
    $basename =~ s/\.clusters$//;
    my $paffile = $outdir.'/'.$basename.'.paf';

    open CLU, '<', $cluster or die "Cant open $cluster: $!\n";
    open PAF, '>', $paffile or die "Can't create $paffile: $!\n";

    while (my $line = <CLU>){

        chomp $line;

        if ($line =~ /^###/){

            ## Cluster data
            my @data = split('; ', $line);
            my $cluster = $data[0];
            my ($qstart,$qend) = split(' to ', $data[1]);
            my ($sstart,$send) = split(' to ', $data[2]);
            my ($len) = $data[3] =~ /size = (\d+)/;

            my $query = $annotations{$qstart}{'contig'};
            my $qm_start = $annotations{$qstart}{'start'};
            my $qm_end = $annotations{$qend}{'end'};
            my $qlen = $lengths{$query};
            my $qstrand = $annotations{$qstart}{'strand'};

            my $subject = $annotations{$sstart}{'contig'};
            my $sm_start = $annotations{$sstart}{'start'};
            my $sm_end = $annotations{$send}{'end'};
            my $slen = $lengths{$subject};
            my $sstrand = $annotations{$sstart}{'strand'};

            my $real_qstart = $qm_start;
            my $real_qend = $qm_end;
            my $real_sstart = $sm_start;
            my $real_ssend = $sm_end;

            ## Adjusting positions and relative strandedness if start > end
            my $rel_strand = '+';
            if ($qm_start > $qm_end){
                $real_qstart = $qm_end;
                $real_qend = $qm_start;
                $rel_strand = '-';
            }

            if ($sm_start > $sm_end){
                $real_sstart = $sm_end;
                $real_ssend = $sm_start;
                $rel_strand = '-';
            }

            ## PAF
            print PAF $subject."\t";
            print PAF $slen."\t";
            print PAF $real_sstart."\t";
            print PAF $real_ssend."\t";
            print PAF $rel_strand."\t";
            print PAF $query."\t";
            print PAF $qlen."\t";
            print PAF $real_qstart."\t";
            print PAF $real_qend."\t";
            print PAF $len."\t";
            print PAF $len."\t";
            print PAF $map_quality."\n";

        }
    }

}