#!/usr/bin/perl
# Pombert lab, 2020

my $name = 'list_maker.pl';
my $version = '0.4.2';
my $updated = '2023-02-03';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;
use File::Path qw(make_path);

my $usage = <<"OPTIONS";
NAME		$name
VERSION		$version
UPDATED		$updated
SYNOPSIS	Creates a list containing all protein-coding genes with their genomic positions from GBF
		(GenBank), GFF, or EMBL files. Also creates .faa files containing protein sequences.

USAGE		$name 
		  -f gbf \\
		  -i *.gbf \\
		  -o LISTS
OPTIONS

$usage .= <<'REGEX';
OPTIONS:
-f (--format)	Input file format: gbf, gff or embl [Default: gbf]
-i (--input)	Input file(s)
-o (--outdir)	Output directory [Default: LISTS]
-r (--regex)	For GFF files; regular expression to parse;  default:

		'^(\S+).*CDS\s+(\d+)\s+(\d+)\s+\.\s+([+-])\s+[012.].*(?:NCBI_GP|Genbank):(\w+\.\d+)'

		# $1 => contig, $2 => start , $3 => end, $4 => strand, $5 => protein
REGEX

die "$usage\n" unless @ARGV;

my %filetypes = (
	'.gbf' => 'gbf',
	'.gbff' => 'gbf',
	'.gb' => '.gbf',
	'.gff' => 'gff',
	'.embl' => 'embl');

my @input_files;
my $format = 'gbf';
my $outdir = 'LIST_MAKER';
my $regex = '^(\S+).*CDS\s+(\d+)\s+(\d+)\s+\.\s+([+-])\s+[012.].*(?:NCBI_GP|Genbank):(\w+\.\d+)';

GetOptions(
	'i|input=s@{1,}' => \@input_files,
	'f|format=s' => \$format,
	'o|outdir=s' => \$outdir,
	'r|regex=s' => \$regex
);

my $list_dir = $outdir."/LISTS";
my $prot_dir = $outdir."/PROT_SEQ";
my $annot_dir = $outdir."/ANNOTATIONS";

my @outdirs = ($list_dir,$prot_dir,$annot_dir);

unless (-d $outdir){
	make_path($outdir,{mode=>0755}) or die("Can't create output directory $outdir: $!\n");
}

for my $dir (@outdirs){
	unless (-d $dir){
		make_path($dir,{mode=>0755}) or die("Can't create output directory $dir: $!\n");
	}
}

print "\n";

foreach my $input_file (@input_files){

	my ($file_name,$dir,$ext) = fileparse($input_file,'\..*');

	print "Creating .list file for $file_name\n";

	open OUT, ">", "$list_dir/$file_name.list" or die "Can't write to output file: $list_dir/$file_name.list\n";
	open PROT, ">", "$prot_dir/$file_name.faa" or die "Can't write to output file: $prot_dir/$file_name.faa\n";
	open ANNOT, ">", "$annot_dir/$file_name.annotations" or die "Can't write to output file: $annot_dir/$file_name.annotations\n";

	if ($filetypes{$ext} eq 'gbf'){
		open GBK, "<", "$input_file" or die "Can't open input file $input_file: $!\n";

		my %location_data;
		
		my $contig;

		my $gene;
		my $CDS;
		my $translate;
		my $section;
		
		my $start;
		my $sections;
		my $end;
		my $strand;
		my $locus;

		my $gene_num = 1;

		my $sequence;

		while (my $line = <GBK>){
			
			chomp $line;

			if ($line =~ /^LOCUS\s+(\S+)/ ){
				$contig = $1;
			}


			## Entering the CDS metadata
			if ($line =~ /^\s{5}CDS/){
				$CDS = 1;
			}
			
			## Accessing CDS metadata
			elsif ($CDS){
				## Leaving CDS metadata
				if ($line !~ /^\s{21}/){
					undef $CDS;
				}

				## Gather CDS metadata
				if ($line =~ /^\s{21}\/locus_tag="(.*)"/){
					$locus = $1;
					$locus =~ s/\W/\_/g;
				}
				elsif ($line =~ /^\s{21}\/product="(.*)"/){
					print ANNOT "$locus\t$1\n";
				}
				elsif ($line =~ /^\s{21}\/translation="([a-zA-Z]+)"*/){
					$translate = 1;
					$sequence .= $1;
				}
				elsif ($translate){
					if ($line =~/^\s{21}([a-zA-Z]+)/){
						$sequence .= $1;
					}
					else{

						my ($start,$end,$strand) = @{$location_data{$locus}};

						print OUT $locus."\t";
						print OUT $contig."\t";
						print OUT $start."\t";
						print OUT $end."\t";
						print OUT $strand."\t";
						print OUT $gene_num."\n";

						print PROT ">$locus\t$contig\t$start\t$end\t$strand\n";
						foreach my $line (unpack("(A60)*",$sequence)){
							print PROT "$line\n";
						}

						undef $translate;
						undef $sequence;
						$gene_num ++;
					}
				}
			}
			
			## Entering the gene metadata
			if ($line =~ /^\s{5}gene\s{12}(?:complement)*\(*<*(\d+)\.\.>*(\d+)/){
				$gene = 1;
				$start = $1;
				$end = $2;
				$strand = "+";
				if ($line =~ /complement/){
					$strand = "-";
				}
			}
			## Entering gene metadata when there is gene on the edge of a circular contig
			elsif ($line =~ /^\s{5}gene\s{12}(?:complement\()*(?:order\()*(?:join\()*(.*)\)/){
				$gene = 1;
				my @sections = split(",",$1);
				my @start = split(/\.\./,$sections[0]);
				my @end = split(/\.\./,$sections[-1]);
				$start = $start[0];
				$end = $end[-1];
				$strand = "+";
				if ($line =~ /complement/){
					$strand = "-";
				}
			}
			elsif ($line =~ /^\s{5}gene\s{12}(?:complement\()*(?:order\()*(?:join\()*(.*)/){
				$gene = 1;
				$section = 1;
				$sections = $1;
			}

			## Accessing gene metadata
			elsif ($gene){
				## Leaving gene metadata
				if ($line !~ /^\s{21}/){
					undef $gene;
				}
				elsif ($section && $line =~ /^s{21}(.*)/){
					$sections .= $1
				}
				elsif ($section && $line =~ /^s{21}\//){
					undef $section;
					my @sections = split(",",$1);
					my @start = split(/\.\./,$sections[0]);
					my @end = split(/\.\./,$sections[-1]);
					$start = $start[0];
					$end = $end[-1];
					$strand = "+";
					if ($line =~ /complement/){
						$strand = "-";
					}
				}
				## Gather gene metadata
				elsif ($line =~ /^\s{21}\/locus_tag="(.*)"/){
					$locus = $1;
					$locus =~ s/\W/\_/g;
					@{$location_data{$locus}} = ($start,$end,$strand);
				}
			}

		}
	}
	elsif($filetypes{$ext} eq 'gff'){
		my %contigs; 
		my %genes;
		my $start;
		my $end;
		open GFF, "<", "$input_file" or die "Can't open input file : $input_file\n";
		while (my $line = <GFF>){
			chomp $line;
			if ($line =~ /^\S+\t\S+\s+\bgene\b\t(\d+)\t(\d+)/){
				$start = $1;
				$end = $2;
			}
			if ($line =~ /$regex/){
				my $contig = $1;
				my $strand = $4;
				my $gene =$5;
				if (!exists $genes{$gene}){
					$genes{$gene} = $gene;
					#$start = sprintf("%10d", $start); $end = sprintf("%10d", $end); ## Preventing number ordering SNAFU
					$contigs{$contig}{$start}[0] = $contig;
					$contigs{$contig}{$start}[1] = $start;
					$contigs{$contig}{$start}[2] = $end;
					$contigs{$contig}{$start}[3] = $strand;
					$contigs{$contig}{$start}[4] = $gene;
				}
			}
		}
		foreach my $cg (sort keys %contigs) { ## We must reorder genes in case they were entered out of order in the GFF files
			my $position = 0;
			foreach my $pos (sort keys %{ $contigs{$cg} }) {
				$position ++;
				print OUT "$contigs{$cg}{$pos}[4]\t"; ## Printing gene
				print OUT "$contigs{$cg}{$pos}[0]\t"; ## Printing contig
				print OUT "$contigs{$cg}{$pos}[1]\t"; ## Printing $start
				print OUT "$contigs{$cg}{$pos}[2]\t"; ## Printing $end
				print OUT "$contigs{$cg}{$pos}[3]\t"; ## Printing $strand
				print OUT "$position\n"; ## Printing $strand
			}
		}
	}
	elsif($filetypes{$ext} eq 'embl'){
		open EMBL, "<", "$input_file" or die "Can't open input file : $input_file\n";
		my $contig = $file_name;
		my $gene = 0;
		my $locus_tag;
		my $strand;
		while (my $line = <EMBL>){
			chomp $line;
			if ($line =~ /FT\s+gene/){ $gene++; }
			elsif ($line =~ /FT\s+\/locus_tag=\"(\w+)\"/){ $locus_tag = $1; }
			elsif ($line =~ /FT\s+CDS/){
				if ($line =~ /complement/){ $strand = '-';}
				else { $strand = '+';}
				my ($start, $end) = $line =~ /\(?(\d+).*\.\.(\d+)/;
				print OUT "$locus_tag"."\t";
				print OUT "$contig"."\t";
				print OUT "$start"."\t";
				print OUT "$end"."\t";
				print OUT "$strand"."\t";
				print OUT "$gene"."\n";
			}
		}
	}
	else{
		if ($format eq 'gbf'){
			open GBK, "<", "$input_file" or die "Can't open input file $input_file: $!\n";
			my $contig;
			my $gene = 0;
			my $locus_tag;
			my $protein_id;
			my $strand;
			my $start; 
			my $end;
			while (my $line = <GBK>){
				
				chomp $line;

				if ($line =~ /^LOCUS\s+(\S+)/ ){
					$contig = $1;
				}
				
				if ($line =~ /^\s+gene\s+.*<(\d+).*>(\d+)/){
					$start = $1;
					$end = $2; 
					$gene++;
				}
				elsif ($line =~ /^\s+\/locus_tag=\"(\w+)\"/){
					$locus_tag = $1;
				}
				elsif ($line =~ /^\s+\/protein_id=\"(\S+)\"/){
					$protein_id = $1;
					if($protein_id =~ /\w+\:(\w+)/){
						$protein_id = $1;
					}
					print LINK "$protein_id\t$locus_tag\n";
				}
				elsif ($line =~ /^\s+CDS/){
					if ($line =~ /complement/){
						$strand = '-';
					}
					else {
						$strand = '+';
					}
					print OUT "$locus_tag"."\t";
					print OUT "$contig"."\t";
					print OUT "$start"."\t";
					print OUT "$end"."\t";
					print OUT "$strand"."\t";
					print OUT "$gene"."\n";
				}
			}
		}
		elsif($format eq 'gff'){
			my %contigs; 
			my %genes;
			my $start;
			my $end;
			open GFF, "<", "$input_file" or die "Can't open input file : $input_file\n";
			while (my $line = <GFF>){
				chomp $line;
				if ($line =~ /^\S+\t\S+\s+\bgene\b\t(\d+)\t(\d+)/){
					$start = $1;
					$end = $2;
				}
				if ($line =~ /$regex/){
					my $contig = $1;
					my $strand = $4;
					my $gene =$5;
					if (!exists $genes{$gene}){
						$genes{$gene} = $gene;
						#$start = sprintf("%10d", $start); $end = sprintf("%10d", $end); ## Preventing number ordering SNAFU
						$contigs{$contig}{$start}[0] = $contig;
						$contigs{$contig}{$start}[1] = $start;
						$contigs{$contig}{$start}[2] = $end;
						$contigs{$contig}{$start}[3] = $strand;
						$contigs{$contig}{$start}[4] = $gene;
					}
				}
			}
			foreach my $cg (sort keys %contigs) { ## We must reorder genes in case they were entered out of order in the GFF files
				my $position = 0;
				foreach my $pos (sort keys %{ $contigs{$cg} }) {
					$position ++;
					print OUT "$contigs{$cg}{$pos}[4]\t"; ## Printing gene
					print OUT "$contigs{$cg}{$pos}[0]\t"; ## Printing contig
					print OUT "$contigs{$cg}{$pos}[1]\t"; ## Printing $start
					print OUT "$contigs{$cg}{$pos}[2]\t"; ## Printing $end
					print OUT "$contigs{$cg}{$pos}[3]\t"; ## Printing $strand
					print OUT "$position\n"; ## Printing $strand
				}
			}
		}
		elsif($format eq 'embl'){
			open EMBL, "<", "$input_file" or die "Can't open input file : $input_file\n";
			print "Working on input file $input_file in $format format\n";
			my $contig = $file_name;
			my $gene = 0;
			my $locus_tag;
			my $strand;
			while (my $line = <EMBL>){
				chomp $line;
				if ($line =~ /FT\s+gene/){ $gene++; }
				elsif ($line =~ /FT\s+\/locus_tag=\"(\w+)\"/){ $locus_tag = $1; }
				elsif ($line =~ /FT\s+CDS/){
					if ($line =~ /complement/){ $strand = '-';}
					else { $strand = '+';}
					my ($start, $end) = $line =~ /\(?(\d+).*\.\.(\d+)/;
					print OUT "$locus_tag"."\t";
					print OUT "$contig"."\t";
					print OUT "$start"."\t";
					print OUT "$end"."\t";
					print OUT "$strand"."\t";
					print OUT "$gene"."\n";
				}
			}
		}
		else{
			print "[E]  Unsupproted annotation file provided ($ext).\n";
			print "\t Supported annotation files: .embl, .gb, .gbff, .gff";
		}
	}

	close OUT;
	close LINK;

}
print "\n";