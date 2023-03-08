#!/usr/bin/perl
## Pombert JF,  Illinois Tech 2020
my $version = '1.4';
my $name = 'queryNCBI.pl';
my $updated = '2021-03-13';

use strict; use warnings; use Getopt::Long qw(GetOptions);

## Usage definition
my $usage = <<"OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Retrieve multifasta files from NCBI's FTP
		TAB/CSV-delimited lists can be generated from NCBI's Genome Assembly and Annotation reports, e.g.:
		1) goto http://www.ncbi.nlm.nih.gov/genome/genomes/159?
		2) click on the Download Table link in the upper right corner 
		
EXAMPLE		${name} \\
		  -l genome_list.csv \\
		  -o DATASETS \\
		  -fa \\
		  -gb \\
		  -gff \\
		  -p \\
		  -cds

OPTIONS:
-l (--list)	TAB/CSV-delimited list from NCBI
-o (--outdir)	Output directory [Default: ./]
-fa (--fasta)	Retrieve fasta files
-gb (--genbank)	Retrieve GenBank annotation files (.gbk; if available)
-gff (--gff3)	Retrieve GFF3 annotation files (.gff; if available)
-p (--protein)	Retrieve protein sequences (.faa; if available)
-cds		Retrieve protein coding sequences (.fna; if available)
OPTIONS
die "\n$usage\n" unless @ARGV;

## Defining options
my $list;
my $outdir = './';
my $fasta;
my $gbk;
my $gff;
my $protein;
my $cds;
GetOptions(
	'l|list=s' => \$list,
	'o|outdir=s' => \$outdir,
	'fa|fasta' => \$fasta,
	'gb|genbank' => \$gbk,
	'gff|gff3' => \$gff,
	'p|protein' => \$protein,
	'cds' => \$cds
);

## Creating output directory
unless (-d $outdir){
	mkdir ($outdir,0755) or die "Can't create folder $outdir: $!\n";
}

## Downloading data from NCBI
my $start = localtime(); my $tstart = time;
system "dos2unix $list"; ## making sure that the line breaks are in UNIX format
open IN, "<", "$list" or die "Can't read file $list: $!\n";

## File information
my $url;
my $accession;
my $organism;
my $genus;
my $species;
my $strain;

while (my $line = <IN>){
	chomp $line;
	$line =~ s/\"//g; ## Removing quotes
	if ($line =~ /^#/){next;} ## Discarding comments
	else {
		my @genome = split(/\t|,/, $line);
		$organism = $genome[0]; ## Organism info
		if ($organism =~ /^(\S+)\s+(\S+)/){ ## Getting rid of random strain tags
			$genus = $1;
			$species = $2;
		}
		
		# Strain info
		$strain = $genome[2];
		$strain =~ s/ /_/g; ## Deleting or substituting problematic characters with underscores
		$strain =~ s/\//_/g;
		$strain =~ s/\|/_/g;
		$strain =~ s/\./_/g;
		$strain =~ s/\(//g;
		$strain =~ s/\)//g;
		
		# Accession info
		my $genbankFTP = $genome[$#genome];
		if ($genbankFTP =~ /.*\/(\S+)$/){$accession = $1;}
		
		## Downloading files
		if ($fasta){
			$url = $genbankFTP.'/'.$accession.'_genomic.fna.gz';
			curl('fasta');
		}
		if ($gbk){
			$url = $genbankFTP.'/'.$accession.'_genomic.gbff.gz';
			curl('gbk');
		}
		if ($gff){
			$url = $genbankFTP.'/'.$accession.'_genomic.gff.gz';
			curl('gff');
		}
		if ($protein){
			$url = $genbankFTP.'/'.$accession.'_protein.faa.gz';
			curl('faa');
		}
		if ($cds){
			$url = $genbankFTP.'/'.$accession.'_cds_from_genomic.fna.gz';
			curl('fna');
		}
	}
}
my $end = localtime(); my $time_taken = time - $tstart;
print "\nDownloads started on: $start\n";
print "Downloads completed on: $end\n";
print "Time elapsed: $time_taken seconds\n";

### Subroutine
sub curl{
	my $ext = $_[0];
	my $file = "$genus".'_'."$species".'_'."$strain".".$ext.gz";
	print "\nDownloading file: $url\n";
	system "curl --progress-bar -L $url -o ${outdir}/$file";
	system "gunzip ${outdir}/$file"; ## decompressing file
}
