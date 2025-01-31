#!/usr/bin/env perl

my $name = 'jgi_to_ncbi_gff.pl';
my $version = '0.1a';
my $updated = '2025-01-31';

use strict;
use warnings;
use File::Basename;
use File::Path qw(make_path);
use Getopt::Long qw(GetOptions);

#########################################################################
### Command line options
#########################################################################

my $usage =<<"USAGE";
NAME        ${name}
VERSION     ${version}
UPDATED     ${updated}
SYNOPSIS    Convert JGI gff files to (pseudo)-NCBI GFF3 format

COMMAND     $name \\
              --gff3 jgi.gff \\
              --out ncbi.gff3 

OPTIONS:
-g (--gff3)     JGI GFF files to convert
-o (--out)      Output file
--version       Show script version
USAGE

unless (@ARGV){
    print "\n$usage\n";
    exit(0);
};

my $gff3;
my $outfile = 'J2N';
my $sc_version;
GetOptions(
    'g|gff3=s' => \$gff3,
    'o|out=s' => \$outfile,
    'version' => \$sc_version
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
### Parse JGI GFF file
#########################################################################

open GFF3, '<', $gff3 or die "Can't read $gff3: $!\n";
open OUT, '>', $outfile or die "Can't create $outfile: $!\n";

my %features;

## Parse features
while (my $line = <GFF3>){

    chomp $line;

    if ($line =~ /^#/){
        next;
    }
    elsif ($line =~ /^\s*$/){
        next;
    }
    elsif ($line =~ /(stop|start)_codon/){
        next;
    }
    else{

        my @data = split("\t", $line);
        my ($seqid) = $data[0] =~ /^(\w+)/;
        my $source = $data[1];
        my $type = $data[2];
        my $start = $data[3];
        my $end = $data[4];
        my $score = $data[5];
        my $strand = $data[6];
        my $phase = $data[7];
        my @attributes = split("; ", $data[8]);

        my $feature_id;
        my $transcriptId;
        my $proteinId;
        my $product;

        foreach my $attribute (@attributes){
            if ($attribute =~ /^name \"(.*)\"/){
                $feature_id = $1;
                $feature_id =~ s/[#!.]/_/g;
                $feature_id =~ s/___/_/g;
                $feature_id =~ s/__/_/g;
            }
            if ($attribute =~ /^transcriptId (.*)/){
                $transcriptId = $1; 
            }
            if ($attribute =~ /^proteinId (.*)/){
                $proteinId = $1; 
            }
            if ($attribute =~ /^(product_name|product) \"(.*)\"/){
                $product = $2;
            }
        }

        ## mRNA data
        if ($type =~ /exon/){

            $features{$seqid}{$feature_id}{'mrna_id'} = $transcriptId;
            my $exons = $start."\t".$end;
            push (@{$features{$seqid}{$feature_id}{'mrna_exons'}}, $exons);
            push (@{$features{$seqid}{$feature_id}{'positions'}}, $start);
            push (@{$features{$seqid}{$feature_id}{'positions'}}, $end);
        }
        elsif ($type =~ /CDS/){

            $features{$seqid}{$feature_id}{'cds_id'} = $proteinId;
            my $exons = $start."\t".$end;
            push (@{$features{$seqid}{$feature_id}{'cds_exons'}}, $exons);
            push (@{$features{$seqid}{$feature_id}{'positions'}}, $start);
            push (@{$features{$seqid}{$feature_id}{'positions'}}, $end);

            my $min;
            if ($start < $end){
                $min = $start;
            }
            else{
                $min = $end;
            }

            unless (exists $features{$seqid}{$feature_id}{'init'}){
                $features{$seqid}{$feature_id}{'init'} = $min;
            }
            else {
                if ($min < $features{$seqid}{$feature_id}{'init'}){
                    $features{$seqid}{$feature_id}{'init'} = $min;
                }
            }

            $features{$seqid}{$feature_id}{'source'} = $source;
            $features{$seqid}{$feature_id}{'score'} = $score;
            $features{$seqid}{$feature_id}{'strand'} = $strand;
            push(@{$features{$seqid}{$feature_id}{'phase'}}, $phase);

            if (defined $product){
                $features{$seqid}{$feature_id}{'product'} = $product;
            }
            else{
                $features{$seqid}{$feature_id}{'product'} = 'undefined product';
            }

        }

    }

}

## Print parsed features
foreach my $seqid (sort (keys %features)){

    my %feat_pos;

    foreach my $feature_id (keys %{$features{$seqid}}){
        my $init = $features{$seqid}{$feature_id}{'init'};
        $feat_pos{$init} = $feature_id;
    }

    # Sorting features by positions
    foreach my $pos (sort { $a <=> $b }(keys %feat_pos)){

        my $feature_id = $feat_pos{$pos};

        my @sorted_coor = sort(@{$features{$seqid}{$feature_id}{'positions'}});
        my $gmin = $sorted_coor[0];
        my $gmax = $sorted_coor[-1];

        # Gene feature
        print OUT $seqid."\t";
        print OUT $features{$seqid}{$feature_id}{'source'}."\t";
        print OUT 'gene'."\t";
        print OUT $gmin."\t";
        print OUT $gmin."\t";
        print OUT $features{$seqid}{$feature_id}{'score'}."\t";
        print OUT $features{$seqid}{$feature_id}{'strand'}."\t";
        print OUT '.'."\t";
        print OUT "ID=gene-$feature_id;Name=$feature_id;locus_tag=$feature_id\n";

        # mRNA feature
        print OUT $seqid."\t";
        print OUT $features{$seqid}{$feature_id}{'source'}."\t";
        print OUT 'mRNA'."\t";
        print OUT $gmin."\t";
        print OUT $gmin."\t";
        print OUT $features{$seqid}{$feature_id}{'score'}."\t";
        print OUT $features{$seqid}{$feature_id}{'strand'}."\t";
        print OUT '.'."\t";
        my $mRNA_id = $features{$seqid}{$feature_id}{'mrna_id'};
        print OUT "ID=rna-$mRNA_id;Parent=gene-$feature_id;Name=$mRNA_id;locus_tag=$feature_id\n";

        # Exon feature(s)
        for my $exons (@{$features{$seqid}{$feature_id}{'mrna_exons'}}){

            my (@pos) = split("\t", $exons);

            print OUT $seqid."\t";
            print OUT $features{$seqid}{$feature_id}{'source'}."\t";
            print OUT 'exon'."\t";
            print OUT $pos[0]."\t";
            print OUT $pos[1]."\t";
            print OUT $features{$seqid}{$feature_id}{'score'}."\t";
            print OUT $features{$seqid}{$feature_id}{'strand'}."\t";
            print OUT '.'."\t";
            my $mRNA_id = $features{$seqid}{$feature_id}{'mrna_id'};
            print OUT "ID=exon-$mRNA_id;Parent=rna-$mRNA_id;Name=$mRNA_id;locus_tag=$feature_id;transcript_id=$mRNA_id\n";
        }

        # CDS feature(s)
        my $c_counter = 0;
        for my $cds (@{$features{$seqid}{$feature_id}{'cds_exons'}}){

            my (@pos) = split("\t", $cds);
            my $cphase = $features{$seqid}{$feature_id}{'phase'}[$c_counter];
            $c_counter++;

            print OUT $seqid."\t";
            print OUT $features{$seqid}{$feature_id}{'source'}."\t";
            print OUT 'CDS'."\t";
            print OUT $pos[0]."\t";
            print OUT $pos[1]."\t";
            print OUT $features{$seqid}{$feature_id}{'score'}."\t";
            print OUT $features{$seqid}{$feature_id}{'strand'}."\t";
            print OUT $cphase ."\t";
            my $cds_id = $features{$seqid}{$feature_id}{'cds_id'};
            my $prod = $features{$seqid}{$feature_id}{'product'};
            print OUT "ID=cds-$cds_id;Parent=rna-$mRNA_id;Name=$cds_id;locus_tag=$feature_id;product=$prod;protein_id=$cds_id\n";
        }

    }
}