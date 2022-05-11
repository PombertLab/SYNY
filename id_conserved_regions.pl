#!/usr/bin/perl
# Pombert Lab 2022

use strict; use warnings; use Getopt::Long qw(GetOptions);

my $name = "id_conserved_regions.pl";
my $version = '0.5a';
my $updated = '2022-05-10';

my $usage = <<"EXIT";
NAME		$name
VERSION		$version
UPDATED		$updated
SYNOPSIS	This scirpt identifies conserved proteins and chromosomal regions between provided species.

USAGE		$name \\
		-l SYNY/LISTS/ \\
		-b SYNY/DIAMOND/ \\
		-o SYNY/SYNTENY/CONSERVED

OPTIONS
-l (--lists)	Directory containing .lists files generated by list_maker.pl
-b (--blasts)	Directory containing DIAMOND blast results
-o (--outdir)	Output directory

EXIT

die("\n$usage\n") unless(@ARGV);

my $lists;
my $blasts;
my $outdir = "CONSERVED";

GetOptions(
	"l|lists=s" => \$lists,
	"b|blasts=s" => \$blasts,
	"o|outdir=s" => \$outdir,
);

unless(-d $outdir){
	mkdir($outdir,0755) or die("Unable to create directory $outdir: $!\n");
}

opendir(LIST,$lists);

my %organisms;
my %protein_bps;
my $org_counter = 0;

foreach my $file (readdir(LIST)){
	if ($file =~ /\.list$/){

		## Obtain the start/end base for each protein as well as each organisms name
		my ($org) = $file =~ /(\w+)\.list/;
		$organisms{$org} = $org_counter;

		open IN, "<", "$lists/$file" or die("Unable to open file $lists/$file: $!\n");
		while(my $line = <IN>){
			chomp($line);
			my ($locus,$accession,$start,$end,undef,undef) = split("\t",$line);
			@{$protein_bps{$org}{$locus}} = ($accession,$start,$end);
			unless($start||$end){
				print($start."\t".$end."\n".$line."\n");
			}
		}
		close IN;
		$org_counter ++;
	}
}

my %blast_hits;
opendir(BLAST,$blasts);
foreach my $file (readdir(BLAST)){
	my %used;
	unless(-d "$blasts/$file"){
		open FILE, "<", "$blasts/$file";
		print "$blasts/$file\n";
		my ($sub_org,$query_org) = $file =~ /^(\w+)_vs_(\w+)\.diamond\.6/;
		while (my $line = <FILE>){
			chomp($line);
			my @data = split("\t",$line);
			my $sub = $data[0];
			my $query = $data[1];
			my $eval = $data[10];

			if ($used{$query}){
				if ($used{$query}[1] > $eval){
					undef($blast_hits{$sub_org}{$query_org}{$used{$query}[0]});
					$blast_hits{$sub_org}{$query_org}{$sub} = $query;
					@{$used{$query}} = ($sub,$eval);
				}
			}
			else{
				@{$used{$query}} = ($sub,$eval);
				$blast_hits{$sub_org}{$query_org}{$sub} = $query;
			}

			if($blast_hits{$sub_org}{$query_org}{$query}){
				if(@{$blast_hits{$sub_org}{$query_org}{$query}}[1] > $eval ){
					# print "Replacing ".@{$blast_hits{$sub_org}{$query_org}{$query}}[0]." with $sub\n";
					@{$blast_hits{$sub_org}{$query_org}{$query}} = ($sub,$eval);
				}
			}
			else{
				@{$blast_hits{$sub_org}{$query_org}{$query}} = ($sub,$eval);
			}
		}
		close FILE;
	}
}

### This chunk of code identifies if a protein is conserved by doing bidirectional analysis of all BLAST results
###
### Example:
###
###	Protein			BLAST result
### Cuniculi_1		(Hc,Ic)
### Hellem_1		(Ch,Ih)
### Intestinalis_1	(Ci,Hi)
###
### If the protein is conserved, then Ch = Ci = Cunculi_1 = Hc = Hi = Hellem_1 = Ic = Ih = Intestinalis_1. To determine
### this, a bidirectional analysis is first performed on Ch and Hc such that if Hc is Hellem_1 and Ch is Cuniculi_1,
### then Hc = Hellem_1 = Ch = Cuniculi_1; a similar analysis is done for Ic and Ci. Next, a bidirectional analysis is
### performed on Ih and Hi such that if Ih is Intestinalis_1 and Hi is Hellem_1, then Ih = Intestinalis_1 = Hi = Hellem_1.
### Combining the two comparisons, the equality Ch = Ci = Cunculi_1 = Hc = Hi = Hellem_1 = Ic = Ih = Intestinalis_1 forms.

my %conserved;
my %unique;

## Iterate over each organism (subject organism)
foreach my $org_sub (sort(keys(%organisms))){

	## File header logistics
	open OUT, ">", "$outdir/$org_sub.conserved_summary" or die "Unable to write to $outdir/$org_sub.conserved_summary: $!\n";

	print OUT "### $org_sub\tHomology Matrix";

	foreach my $org_name (sort(keys(%organisms))){
		unless ($org_name eq $org_sub){
			print OUT "\t".$org_name;
		}
	}

	print OUT "\n";

	my $previous_accession;
	foreach my $sub_prot (sort(keys(%{$protein_bps{$org_sub}}))){
		my $binary = "";
		foreach my $org_query (sort(keys(%organisms))){
			unless($org_sub eq $org_query){
				if($blast_hits{$org_sub}{$org_query}{$sub_prot}){
					if($blast_hits{$org_query}{$org_sub}{$sub_prot}){
						if($blast_hits{$org_query}{$org_sub}{$sub_prot}[0] eq $blast_hits{$org_sub}{$org_query}{$sub_prot}){
							$binary .= "2";
						}
						else{
							$binary .= "1";
						}
					}
					else{
						$binary .= "1";
					}
				}
				else{
					$binary .= "0";
				}
			}
			else{
				$binary .= "-";
			}
		}

		my $pos = 0;
		foreach my $org_query (sort(keys(%organisms))){
			unless($org_query eq $org_sub){
				$binary .= " ";
				if($blast_hits{$org_sub}{$org_query}{$sub_prot}){
					my $query_prot = $blast_hits{$org_sub}{$org_query}{$sub_prot};
					foreach my $org_comp (sort(keys(%organisms))){
						unless($org_comp eq $org_query){
							if($blast_hits{$org_query}{$org_comp}{$query_prot}){
								if($blast_hits{$org_comp}{$org_query}{$query_prot}){
									if($blast_hits{$org_comp}{$org_query}{$query_prot}[0] eq $blast_hits{$org_query}{$org_comp}{$query_prot}){
										$binary .= "2";
									}
									else{
										$binary .= "1";
									}
								}
								else{
									$binary .= "1";
								}
							}
							else{
								$binary .= "0";
							}
						}
						else{
							$binary .= "-";
						}
					}
				}
				else{
					$binary .= ("0"x$pos)."-".("0"x(scalar(%organisms)-$pos-1));
				}
			}
			$pos ++;
		}

		if($previous_accession){
			if ($protein_bps{$org_sub}{$sub_prot}[0] ne $previous_accession){
				print OUT "## ".$protein_bps{$org_sub}{$sub_prot}[0]."\n";
			}
		}
		else{
			print OUT "## ".$protein_bps{$org_sub}{$sub_prot}[0]."\n";
		}

		$previous_accession = $protein_bps{$org_sub}{$sub_prot}[0];

		print OUT $sub_prot."\t".$binary;

		if ($binary !~ /[01]/){
			$conserved{$org_sub}{$sub_prot} = 0;
		}

		if ($binary !~ /[12]/){
			$unique{$org_sub}{$sub_prot} = 0;
		}

		foreach my $org_query (sort(keys(%organisms))){
			unless($org_sub eq $org_query){
				if($blast_hits{$org_sub}{$org_query}{$sub_prot}){
					print OUT "\t".$blast_hits{$org_sub}{$org_query}{$sub_prot};
				}
				else{
					print OUT "\t"."-";
				}
			}
		}
		print OUT "\n";
	}

	close OUT;

}

foreach my $org (sort(keys(%protein_bps))){

	open OUT, ">", "$outdir/$org.conserved" or die "Unable to write to $outdir/$org.conserved: $!\n";

	my $previous_accession;
	my $previous_protein;

	print OUT "#"x100;
	print OUT "\n## Contig Level\n";
	print OUT "#"x100;

	print OUT "\n## Contig Accession\tFirst Conserved Protein Location\tLast Conserved Protein Location\n";

	foreach my $protein (sort(keys(%{$conserved{$org}}))){

		unless($previous_accession){
			$previous_accession = $protein_bps{$org}{$protein}[0];
			print OUT "\n".$previous_accession."\t".$protein_bps{$org}{$protein}[1]."-".$protein_bps{$org}{$protein}[2];
		}
		elsif($previous_accession ne $protein_bps{$org}{$protein}[0]){
			print OUT "\t".$protein_bps{$org}{$previous_protein}[2]."-".$protein_bps{$org}{$previous_protein}[2]."\n";
			$previous_accession = $protein_bps{$org}{$protein}[0];
			print OUT $previous_accession."\t".$protein_bps{$org}{$protein}[1]."-".$protein_bps{$org}{$protein}[2];
		}
		$previous_protein = $protein;

	}

	print OUT "\t".$protein_bps{$org}{$previous_protein}[1]."-".$protein_bps{$org}{$previous_protein}[2]."\n\n";

	print OUT "#"x100;
	print OUT "\n## Protein Level\n";
	print OUT "#"x100;
	print OUT "\n##Contig Accession\tProtein Location\t".$org. "locus";
	foreach my $org_query (sort(keys(%{$blast_hits{$org}}))){
		unless($org eq $org_query){
			print OUT "\t".$org_query." locus";
		}
	}
	print OUT "\n\n";

	foreach my $protein (sort(keys(%{$conserved{$org}}))){
		
		print OUT $protein_bps{$org}{$protein}[0]."\t".$protein_bps{$org}{$protein}[1]."-".$protein_bps{$org}{$protein}[2]."\t".$protein;

		foreach my $org_query (sort(keys(%{$blast_hits{$org}}))){
			unless($org eq $org_query){
				print OUT "\t".$blast_hits{$org}{$org_query}{$protein};
			}
		}
		print OUT "\n"

	}

}

foreach my $org (sort(keys(%protein_bps))){
	open OUT, ">", "$outdir/$org.unique" or die "Unable to write to $outdir/$org.unique: $!\n";

	print OUT "#"x100;
	print OUT "\n## Protein Level\n";
	print OUT "#"x100;
	print OUT "\n##Contig Accession\tProtein Location\t".$org. "locus\n\n";

	foreach my $protein (sort(keys(%{$unique{$org}}))){
		
		print OUT $protein_bps{$org}{$protein}[0]."\t".$protein_bps{$org}{$protein}[1]."-".$protein_bps{$org}{$protein}[2]."\t".$protein."\n";

	}

	close OUT;
}