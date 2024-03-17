#!/usr/bin/env perl
# Pombert lab, 2020

my $name = 'list_maker.pl';
my $version = '0.5.3a';
my $updated = '2024-03-17';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;
use File::Path qw(make_path);
use PerlIO::gzip;

my $usage = <<"OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Creates a list containing all protein-coding genes with their genomic positions from GBF
		GenBank GBF/GBFF annotation files. Also creates .faa files containing protein sequences.

USAGE		${name}
		  -i *.gbf \\
		  -o LISTS
OPTIONS

$usage .= <<'REGEX';
OPTIONS:
-i (--input)	Input file(s) (GZIP supported; File type determined by file extension)
-o (--outdir)	Output directory [Default: LIST_MAKER]
REGEX

die "$usage\n" unless @ARGV;

my %filetypes = (
	'gbf' => 'gbf',
	'gbff' => 'gbf',
	'gb' => '.gbf',
	'gff' => 'gff',
	'embl' => 'embl');

my @input_files;
my $outdir = 'LIST_MAKER';

GetOptions(
	'i|input=s@{1,}' => \@input_files,
	'o|outdir=s' => \$outdir
);


###################################################################################################
## Output dir/subdir creation and setup
###################################################################################################

my $list_dir = $outdir."/LISTS";
my $prot_dir = $outdir."/PROT_SEQ";
my $genome_dir = $outdir."/GENOME";

my @outdirs = ($list_dir,$prot_dir,$genome_dir);

unless (-d $outdir){
	make_path($outdir,{mode=>0755}) or die("Can't create output directory $outdir: $!\n");
}

for my $dir (@outdirs){
	unless (-d $dir){
		make_path($dir,{mode=>0755}) or die("Can't create output directory $dir: $!\n");
	}
}

###################################################################################################
## Parsing input files
###################################################################################################

print "\n";

foreach my $input_file (@input_files){

	my $file_name = basename($input_file);
	my @file_data = split('\.',$file_name);
	my $file_prefix = $file_data[0];
	my $diamond = "<";
	my $ext = $file_data[-1];

	if ($file_data[-1] eq 'gz'){ 
		$diamond = "<:gzip";
		$ext = $file_data[-2];
	}

	print "Creating .list file for $file_name\n";


	my $annotation_list = "$list_dir/$file_prefix.list";
	my $fasta_proteins = "$prot_dir/$file_prefix.faa";
	my $fasta_genome = "$genome_dir/$file_prefix.fasta";

	open OUT, ">", $annotation_list or die "Can't create $annotation_list: $!\n";
	open PROT, ">", $fasta_proteins or die "Can't create $fasta_proteins: $!\n";
	open GEN, ">", $fasta_genome or die "Can't create $fasta_genome: $!\n";

	## GenBank GBF/GBFF format
	if ($filetypes{$ext} eq 'gbf'){

		open GBK, $diamond, $input_file or die "Can't open input file $input_file: $!\n";

		my %location_data;
		my %features;
		my %genome;

		my $seq_flag;
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

		my $incomplete_product_name;

		while (my $line = <GBK>){
			
			chomp $line;

			if ($line =~ /^LOCUS\s+(\S+)/ ){
				$contig = $1;
				## if the contig is not named with a unique LOCUS
				## add the file_prefix in front to prevent naming clashes 
				if ($contig =~ /chromosome|contig/i){
					$contig = $file_prefix.'_'.$contig;
				}
			}

			if ($line =~ /^ORIGIN/){
				$seq_flag = 1;
			}
			elsif ($line =~ /^\/\//){
				$seq_flag = undef;
			}

			if ($seq_flag){
				unless ($line =~ /^ORIGIN/){
					$line =~ s/\s//g;
					$line =~ s/\d//g;
					$genome{$contig} .= uc($line);
				}
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

				## Checking for product names, including accross multiple lines
				elsif ($line =~ /^\s{21}\/product="(.*)"/){
					$features{$locus}{'product'} = $1;
				}
				elsif ($line =~ /^\s{21}\/product="(.*)/){
					$features{$locus}{'product'} = $1;
					$incomplete_product_name = $1;
				}
				elsif ((defined $incomplete_product_name) && ($line =~ /^\s{21}(.*)\"$/)){
					$features{$locus}{'product'} .= " $1";
					$incomplete_product_name = undef;
				}
				elsif ((defined $incomplete_product_name) && ($line =~ /^\s{21}(.*)/)){
					$features{$locus}{'product'} .= " $1";
					$incomplete_product_name .= " $1";
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
						print OUT $gene_num."\t";

						if (!defined $features{$locus}{'product'}){
							$features{$locus}{'product'} = 'undefined product in accession';
						}
						print OUT $features{$locus}{'product'}."\n";

						print PROT ">$locus \[$features{$locus}{'product'}\]\n"; #\t$contig\t$start\t$end\t$strand\n";
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

		### Creating genome FASTA from GBFF
		for my $sequence (sort(keys %genome)){
			print GEN ">$sequence\n";
			my @data = unpack ("(A60)*", $genome{$sequence});
			while (my $tmp = shift@data){
				print GEN "$tmp\n";
			}
		}

	}

	### Other formats
	else {
		print "Unrecognized file type: $filetypes{$ext}\n";
		exit;
	}

	# GFF
	# elsif ($filetypes{$ext} eq 'gff'){
	# 	my %contigs; 
	# 	my %genes;
	# 	my $start;
	# 	my $end;
	# 	open GFF, $diamond, "$input_file" or die "Can't open input file : $input_file\n";
	# 	while (my $line = <GFF>){
	# 		chomp $line;
	# 		if ($line =~ /^\S+\t\S+\s+\bgene\b\t(\d+)\t(\d+)/){
	# 			$start = $1;
	# 			$end = $2;
	# 		}
	# 		if ($line =~ /$regex/){
	# 			my $contig = $1;
	# 			my $strand = $4;
	# 			my $gene = $5;
	# 			my $product = $6;
	# 			if (!exists $genes{$gene}){
	# 				$genes{$gene} = $gene;
	# 				#$start = sprintf("%10d", $start); $end = sprintf("%10d", $end); ## Preventing number ordering SNAFU
	# 				$contigs{$contig}{$start}[0] = $contig;
	# 				$contigs{$contig}{$start}[1] = $start;
	# 				$contigs{$contig}{$start}[2] = $end;
	# 				$contigs{$contig}{$start}[3] = $strand;
	# 				$contigs{$contig}{$start}[4] = $gene;
	# 				$contigs{$contig}{$start}[5] = $product;
	# 			}
	# 		}
	# 	}
	# 	foreach my $cg (sort keys %contigs) { ## We must reorder genes in case they were entered out of order in the GFF files
	# 		my $position = 0;
	# 		foreach my $pos (sort keys %{ $contigs{$cg} }) {
	# 			$position ++;
	# 			print OUT "$contigs{$cg}{$pos}[4]\t"; ## Printing gene
	# 			print OUT "$contigs{$cg}{$pos}[0]\t"; ## Printing contig
	# 			print OUT "$contigs{$cg}{$pos}[1]\t"; ## Printing $start
	# 			print OUT "$contigs{$cg}{$pos}[2]\t"; ## Printing $end
	# 			print OUT "$contigs{$cg}{$pos}[3]\t"; ## Printing $strand
	# 			print OUT "$position\t";
	# 			print OUT "$contigs{$cg}{$pos}[5]\n"; ## Printing $product
	# 		}
	# 	}
	# }
	# # EMBL
	# elsif ($filetypes{$ext} eq 'embl'){
	# 	open EMBL, $diamond, "$input_file" or die "Can't open input file : $input_file\n";
	# 	my $contig = $file_name;
	# 	my $gene = 0;
	# 	my $locus_tag;
	# 	my $product;
	# 	my $strand;
	# 	while (my $line = <EMBL>){
	# 		chomp $line;
	# 		if ($line =~ /^ID\s+(\w+);/){ $contig = $1 }
	# 		elsif ($line =~ /FT\s+gene/){ $gene++; }
	# 		elsif ($line =~ /FT\s+\/locus_tag=\"(\w+)\"/){ $locus_tag = $1; }
	# 		elsif ($line =~ /FT\s+\/product=\"(.*)\"/){ $product = $1; }
	# 		elsif ($line =~ /FT\s+CDS/){
	# 			if ($line =~ /complement/){ $strand = '-';}
	# 			else { $strand = '+';}
	# 			my ($start, $end) = $line =~ /\(?(\d+).*\.\.(\d+)/;
	# 			print OUT $locus_tag."\t";
	# 			print OUT $contig."\t";
	# 			print OUT $start."\t";
	# 			print OUT $end."\t";
	# 			print OUT $strand."\t";
	# 			print OUT $gene."\t";
	# 			print OUT $product."\n";
	# 		}
	# 	}
	# }

	# close OUT;

}
print "\n";