#!/usr/bin/env perl
# Pombert lab, 2020

my $name = 'list_maker.pl';
my $version = '0.5.8b';
my $updated = '2024-05-23';

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
SYNOPSIS	Creates a list containing all protein-coding genes with their genomic positions from
		GenBank GBF/GBFF annotation files. Also creates .faa files containing protein sequences.

USAGE		${name}
		  -i *.gbff \\
		  -o LISTS \\
		  -x CATOPV '^CAJUZD'
		  -verbose
OPTIONS

$usage .= <<'REGEX';
OPTIONS:
-i (--input)	Input file(s) (GZIP supported; File type determined by file extension)
-o (--outdir)	Output directory [Default: LIST_MAKER]
-m (--minsize)	Keep only contigs larger or equal to specified size (in bp) [Default: 1]
-n (--include)	Select contigs with names from input text file(s) (one name per line); i.e. excludes everything else
-r (--ranges)	Select contigs with subranges from input text file(s): name start end
-x (--exclude)	Exclude contigs with names matching the provided regular expression(s)
-v (--verbose)	Add verbosity
REGEX
die "$usage\n" unless @ARGV;

my %filetypes = (
	'gbf' => 'gbf',
	'gbff' => 'gbf',
	'gb' => '.gbf',
	'gbk' => '.gbf'
);

my @input_files;
my $outdir = 'LIST_MAKER';
my $minsize = 1;
my @included;
my @ranges;
my @regexes;
my $verbose;
GetOptions(
	'i|input=s@{1,}' => \@input_files,
	'o|outdir=s' => \$outdir,
	'm|minsize=i' => \$minsize,
	'n|included=s@{0,}' => \@included,
	'r|ranges=s@{0,}' => \@ranges,
	'x|exclude=s@{0,}' => \@regexes,
	'v|verbose' => \$verbose
);

###################################################################################################
## Output dir/subdir creation and setup
###################################################################################################

my $list_dir = $outdir.'/LISTS';
my $seq_dir = $outdir.'/SEQUENCES';
my $prot_dir = $seq_dir.'/PROTEINS';
my $genome_dir = $seq_dir.'/GENOMES';

my @outdirs = ($list_dir,$seq_dir,$prot_dir,$genome_dir);

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

my %included;
if (@included){
	my @tmp = @included;
	while (my $inc = shift @tmp){
		open INC, '<', $inc or die "Can't open $inc: $!\n";
		while (my $line = <INC>){
			chomp $line;
			if ($line =~ /^#/){
				next;
			}
			elsif ($line =~ /^(\w+)/){
				my $contig_to_keep = $1;
				$included{$contig_to_keep} = 1;
			}
		}
		close INC;
	}
}

my %ranges;
if (@ranges){
	my @tmp = @ranges;
	while (my $inc = shift @tmp){
		open RNG, '<', $inc or die "Can't open $inc: $!\n";
		while (my $line = <RNG>){
			chomp $line;
			$line =~ s/^\s+//; ## Getting rid of spaces (if any) as the start of the line
			if ($line =~ /^#/){
				next;
			}
			elsif ($line eq ''){
				next;
			}
			else {
				my @data = split(/\s+/, $line);
				my ($cname) = $data[0] =~ /^(\w+)/;
				my $cstart = $data[1];
				my $cend = $data[2];
				push (@{$ranges{$cname}}, "$cstart;$cend");
			}
		}
		close RNG;
	}
}

foreach my $input_file (@input_files){

	my $file_name = basename($input_file);
	my @file_data = split('\.',$file_name);
	my ($file_prefix) = $file_data[0] =~ /^(\w+)/;
	my $diamond = "<";
	my $ext = $file_data[-1];

	if ($file_data[-1] eq 'gz'){ 
		$diamond = "<:gzip";
		$ext = $file_data[-2];
	}

	if ($verbose){
		print "\nCreating .list file for $file_name\n";
	}

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
		my %genome;
		my %isoform; ## Keeping track of isoforms with identical locus tags

		my $seq_flag;
		my $contig;
		my $contig_size;
		my $contig_sum = 0;
		my $contig_kept_counter = 0;
		my $contig_kept_sum = 0;
		my $contig_discarded_counter = 0;
		my $contig_discarded_sum = 0;

		my $gene;
		my $CDS;
		my $product;
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
		my %excluded_seq = ();

		while (my $line = <GBK>){
			
			chomp $line;

			if ($line =~ /^LOCUS\s+(\S+).*?(\d+) bp/ ){

				$contig = $1;
				$contig_size = $2;
				$contig_sum += $contig_size;

				if (@regexes){
					for my $regex (@regexes){
						if ($contig =~ /$regex/){
							$excluded_seq{$contig} = $contig_size;
						}
					}
				}

				if ($verbose){
					if ($contig_size >= $minsize){
						unless (exists $excluded_seq{$contig}){
							$contig_kept_counter++;
							$contig_kept_sum += $contig_size;
							print ' '.$contig."\t".$contig_size."\t". ">= $minsize bp"."\t"."kept\n";
						}
						else {
							$contig_discarded_counter++;
							$contig_discarded_sum += $contig_size;
							print ' '.$contig."\t".$contig_size."\t"."regex match"."\t"."discarded\n";
						}

					}
					else {
						$contig_discarded_counter++;
						$contig_discarded_sum += $contig_size;
						unless (exists $excluded_seq{$contig}){
							print ' '.$contig."\t".$contig_size."\t"."< $minsize bp"."\t"."discarded\n";
						}
						else {
							print ' '.$contig."\t".$contig_size."\t"."< $minsize bp and regex match"."\t"."discarded\n";
						}
					}
				}

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
					$product = $1;
				}
				elsif ($line =~ /^\s{21}\/product="(.*)/){
					$product = $1;
					$incomplete_product_name = $1;
				}
				elsif ((defined $incomplete_product_name) && ($line =~ /^\s{21}(.*)\"$/)){
					$product .= " $1";
					$incomplete_product_name = undef;
				}
				elsif ((defined $incomplete_product_name) && ($line =~ /^\s{21}(.*)/)){
					$product .= " $1";
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

						my $start,
						my $end;
						my $strand;
						if (exists $location_data{$locus}){
							($start,$end,$strand) = @{$location_data{$locus}};
						}
						else{
							next;
						}

						unless (exists $isoform{$locus}){ ## Keeping only the first isoform

							if (@regexes){
								for my $regex (@regexes){
									if ($contig =~ /$regex/){
										next;
									}
								}
							}

							if (@included){
								if (!exists $included{$contig}){
									next;
								}
							}

							if (@ranges){
								my @subranges = @{$ranges{$contig}};
								my $check_range = 0;
								for my $subr (@subranges){
									my ($sstart,$send) = split(';', $subr);
									if (($start >= $sstart) && ($end <= $send)){
										$check_range = 1;
									}
								}
								if ($check_range == 0){
									next;
								}
							}

							if ($contig_size >= $minsize){

								unless (exists $excluded_seq{$contig}){

									$isoform{$locus} = '';
									print OUT $locus."\t";
									print OUT $contig."\t";
									print OUT $start."\t";
									print OUT $end."\t";
									print OUT $strand."\t";
									print OUT $gene_num."\t";

									$product = "undefined product in accession" unless defined $product;
									print OUT $product."\n";

									print PROT ">$locus \[$product\]\n"; #\t$contig\t$start\t$end\t$strand\n";
									foreach my $line (unpack("(A60)*",$sequence)){
										print PROT "$line\n";
									}

								}

							}

							undef $translate;
							undef $sequence;
							undef $product;
							$gene_num ++;

						}

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
				## Gather gene metadata from locus tag
				elsif ($line =~ /^\s{21}\/locus_tag="(.*)"/){
					$locus = $1;
					$locus =~ s/\W/\_/g;
					@{$location_data{$locus}} = ($start,$end,$strand);
				}
				## or from gene_id (if no locus_tag)
				elsif ($line =~ /^\s{21}\/db_xref="GeneID:(\d+)"/){
					$locus = $1;
					@{$location_data{$locus}} = ($start,$end,$strand);
				}
			}

		}

		if ($verbose){
			my $total_contigs = $contig_kept_counter + $contig_discarded_counter;
			print "\n".'Total contigs: '."\t\t\t".$total_contigs."\t".$contig_sum." bp\n";
			print "Contigs kept (>= $minsize bp): "."\t".$contig_kept_counter."\t".$contig_kept_sum." bp\n";
			print "Contigs discarded (< $minsize bp and/or regex): "."\t".$contig_discarded_counter."\t".$contig_discarded_sum." bp\n";
		}

		### Creating genome FASTA from GBFF
		for my $sequence (sort(keys %genome)){

			my $contig_len = length($genome{$sequence});

			## Keeping only contig larger than minsize
			if ($contig_len >= $minsize){

				## Keeping only contigs from provided list
				if (%included){
					if (!exists $included{$sequence}){
						next;
					}
				}

				if (%ranges){
					if (!exists $ranges{$sequence}){
						next;
					}
				}

				## Excluding sequences matching requested regexes
				unless (exists $excluded_seq{$sequence}){

					my $genome_seq = $genome{$sequence};

					my $subrange;
					if (%ranges){
						if (exists $ranges{$sequence}){
			
							my @subranges = @{$ranges{$sequence}};
							my $rng_counter = 0;
							my $rng_total = scalar (@subranges);

							for my $subr (@subranges){

								$rng_counter++;

								my ($cstart,$cend) = split(';', $subr);
								my $sublen = $cend - $cstart + 1;

								$genome_seq = substr($genome{$sequence}, $cstart - 1, $sublen);
								$subrange = 'subrange: ';
								$subrange .= $cstart;
								$subrange .= ' - ';
								$subrange .=  $cend;

								if ($rng_total == 1){
									print GEN ">$sequence \[$subrange\]\n";
								}
								else {
									print GEN ">${sequence}_${rng_counter} \[${subrange}\]\n";
								}

								my @data = unpack ("(A60)*", $genome_seq);
								while (my $tmp = shift@data){
									print GEN "$tmp\n";
								}

							}

						}
					}

					else {
						print GEN ">$sequence\n";
						my @data = unpack ("(A60)*", $genome_seq);
						while (my $tmp = shift@data){
							print GEN "$tmp\n";
						}
					}

				}
			}

		}

	}

	### Other formats
	else {
		if ($verbose){
			print "Unrecognized file type: $filetypes{$ext}\n";
		}
		exit;
	}
}

if ($verbose){
	print "\n";
}