#!/usr/bin/env perl
## Pombert Lab, 2025

my $name = 'paf_minsize.pl';
my $version = '0.1a';
my $updated = '2025-03-26';

use strict;
use warnings;
use File::Basename;
use Cwd qw(abs_path);
use Getopt::Long qw(GetOptions);

#########################################################################
### Command line options
#########################################################################

my $usage =<<"USAGE";
NAME        ${name}
VERSION     ${version}
UPDATED     ${updated}
SYNOPSIS    Filters out alignments smaller than size n from PAF files

EXAMPLE     ${name} \\
              -p *.paf \\
              -m 5000 \\
              -o PAF

PREREQS     Minimap2:     https://github.com/lh3/minimap2

OPTIONS:
-p (--paf)      PAF files to parse
-m (--minsize)  Minimum alignment size to keep [Default: 1000]
-o (--outdir)   Output directory [Default: ./PAF]
-v (--version)  Show script version
USAGE

unless (@ARGV){
	print "\n$usage\n";
	exit(0);
};

my @commands = @ARGV;

my @paf;
my $minsize = 1000;
my $outdir = './PAF_MINSIZE';
my $sc_version;
GetOptions(
    'p|paf=s@{1,}' => \@paf,
    'm|minsize=i' => \$minsize,
    'o|outdir=s' => \$outdir,
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
### Output dir/subdirs
#########################################################################

$outdir = abs_path($outdir);
unless (-d $outdir){
    mkdir ($outdir, 0755) or die "Can't create $outdir: $!\n";
}

#########################################################################
### Parsing PAF files
#########################################################################

foreach my $paf (@paf){

    my $basename = fileparse($paf);
    $basename =~ s/\.paf$//;
    my $outfile = $outdir.'/'.$basename.'.m'.$minsize.'.paf';
    print $outfile."\n";

    open IN, '<', $paf or die "Can't open $paf: $!\n";
    open OUT, '>', $outfile or die "Can't create $outfile: $!\n";

    while (my $line = <IN>){

        chomp $line;

        # Skip comments
        if ($line =~ /^#/){
            next;
        }

        # Skip blank lines
        elsif ($line =~ /^$/){
            next;
        }

        else {

            my @data = split(/\t/, $line);

            # 1	string	Query sequence name
            # 2	int	Query sequence length
            # 3	int	Query start coordinate (0-based)
            # 4	int	Query end coordinate (0-based)

            my $q_start = $data[2];
            my $q_end = $data[3];

            my $a_len = $q_end - $q_start + 1;

            if ($a_len >= $minsize){
                print OUT $line."\n";
            }
        }

    } 

    close IN;
    close OUT;

}