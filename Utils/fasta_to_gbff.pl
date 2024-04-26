#!/usr/bin/env perl
## Pombert Lab, 2024

my $name = 'fasta_to_gbff.pl';
my $version = '0.1';
my $updated = '2024-04-26';

use strict;
use warnings;
use File::Basename;
use File::Path qw(make_path);
use PerlIO::gzip;
use Getopt::Long qw(GetOptions);

my $usage =<<"USAGE";
NAME        ${name}
VERSION     ${version}
UPDATED     ${updated}
SYNOPSIS    Converts FASTA files to GenBank GBFF format (without annotations)

COMMAND     $name \\
              --fasta *.fasta
              --outdir ./
              --gzip

OPTIONS:
-f (--fasta)    FASTA file to convert
-o (--outdir)   Output directory [Default: GBFF]
-g (--gzip)     Compress the GBFF output files
-v (--verbose)  Add verbosity
USAGE
die "\n$usage\n" unless @ARGV;

my @fasta;
my $outdir = 'GBFF';
my $gzip_flag;
my $verbose;
GetOptions(
    'f|fasta=s@{1,}' => \@fasta,
    'o|outdir=s' => \$outdir,
    'g|gzip' => \$gzip_flag,
    'v|verbose' => \$verbose
);

### Pigz check
my $gzip_tool = 'gzip';
my $pigz_check = `echo \$(command -v pigz)`;
if ($pigz_check =~ /pigz/){
    $gzip_tool = 'pigz';
}

### Outdir
unless (-d $outdir){
    make_path($outdir,{mode=>0755}) or die "Can't create $outdir: $!\n";
}

### Parsing fasta file(s)
while (my $fasta = shift@fasta){

    if ($verbose){
        print 'Parsing           : '.$fasta."\n";
    }

    ## checking for gzip file extension
    my $gzip = '';
    if ($fasta =~ /.gz$/){
        $gzip = ':gzip';
    }

    open FASTA, "<$gzip", $fasta or die "Can't read $fasta: $!\n";
    my ($basename,$path) = fileparse($fasta);
    $basename =~ s/\.gz$//;
    $basename =~ s/\.\w+$//;
    my $outfile = $outdir.'/'.$basename.'.gbff';
    open GBFF, '>', $outfile or die "Can't create $outfile: $!\n";

    my %sequences;
    my $header;
    while (my $line = <FASTA>){

        chomp $line;
        if ($line =~ /^>(\w+)/){
            $header = $1;
        }
        else {
            $sequences{$header} .= $line;
        }

    }

    for my $sequence (sort(keys %sequences)){

        ## GBFF header
        my $length = length $sequences{$sequence};
        print GBFF "LOCUS       $sequence           $length bp    DNA"."\n";
        print GBFF "DEFINITION  $sequence"."\n";
        print GBFF "ACCESSION   $sequence"."\n";
        print GBFF "VERSION     $sequence"."\n";
        print GBFF "KEYWORDS    ."."\n";
        print GBFF "SOURCE      ."."\n";
        print GBFF "ORGANISM    ."."\n";
        print GBFF "            ."."\n";
        print GBFF "FEATURES             Location/Qualifiers"."\n";

        ## GBFF sequence
        print GBFF "ORIGIN      "."\n";
        my @sequence = unpack ("(A60)*", lc($sequences{$sequence}));
        my $seq_counter = 1;

        while (my $wide_60 = shift@sequence){

            my @wide_10 = unpack ("(A10)*", $wide_60);
            my $seq_counter_len = length($seq_counter);
            my $spacer_len = 9 - $seq_counter_len;
            my $spacer = ' ' x $spacer_len;

            print GBFF $spacer.$seq_counter;
            while (my $wide_10 = shift@wide_10){
                print GBFF ' '.$wide_10;
            }
            print GBFF "\n";

            $seq_counter += 60;

        }

        print GBFF '//'."\n";

    }

    ## Closing filehandles
    if ($gzip eq ':gzip'){
        binmode FASTA, ":gzip(none)";
    }

    close FASTA;
    close GBFF;

    ## Compress GFF files
    if ($gzip_flag){
        if ($verbose){
            print "Compressing ($gzip_tool): ".$outfile."\n";
        }
        system ("$gzip_tool $outfile");
    }

}
