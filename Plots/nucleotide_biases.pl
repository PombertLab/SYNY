#!/usr/bin/perl
## Pombert Lab, 2022
my $name = 'nucleotide_biases.pl';
my $version = '0.7g';
my $updated = '2024-05-22';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;
use File::Path qw(make_path);
use Text::Roman qw(:all);

my $usage = <<"OPTIONS";
NAME		$name
VERSION		$version
UPDATED		$updated
SYNOPSIS	Generates tab-delimited sliding windows of GC, AT, purine and pyrimidine
		 distributions for easy plotting with MS Excel or other tool. Can also generate
		 coordinate files for Circos.

COMMAND		$name \\
		  -fasta *.fasta \\
		  -outdir output_directory \\
		  -winsize 1000 \\
		  -step 500 \\
		  -minsize 1 \\
		  -reference CCMP1205 \\
		  -gap 0 \\
		  -custom_preset chloropicon

OPTIONS (Main)
-f (--fasta)		Fasta file(s) to process
-o (--outdir)		Output directory [Default: ntBiases]
-w (--winsize)		Sliding window size [Default: 10000]
-s (--step)		Sliding window step [Default: 5000]
-m (--minsize)	Minimum contig size (in bp) [Default: 1]
-n (--ncheck)		Check for ambiguous/masked (Nn) nucleotides
-t (--tsv)		Output tab-delimited files (e.g. for excel plotting)

OPTIONS (Circos data files options)
-r (--reference)	Genome reference for Circos plotting
-g (--gap)		Default gap links file for Circos plotting [Default: 0]
-u (--unit)		Size unit (Kb or Mb) [Default: Mb]
-l (--labels)	Contig label type: mixed (arabic + roman numbers), arabic, roman, or names [Default: mixed]
-label_size		Contig label size [Default: 36]
-label_font		Contig label font [Default: bold]
-custom_file		Load custom colors from file
-list_preset		List available custom color presets
-custom_preset		Use a custom color preset; e.g.
			# chloropicon - 20 colors - Lemieux et al. (2019) https://pubmed.ncbi.nlm.nih.gov/31492891/
			# encephalitozoon - 11 colors - Pombert et al. (2012) https://pubmed.ncbi.nlm.nih.gov/22802648/
			# presets can be added to the cc_colors subroutine
-noticks		Comment out ticks.conf in Circos configuration file
-max_ticks		Set max number of ticks [Default: 5000]
-max_ideograms		Set max number of ideograms [Default: 200]
-max_links		Set max number of links [Default: 75000]
-max_points_per_track	Set max number of points per track [Default: 75000]
-zdepth			Set links zdepth in ruleset [Default: 50]
-clusters		Color by clusters [Default: off]
-no_biases		Skip nucleotide bias subplots in Circos configuration files
OPTIONS
die "\n$usage\n" unless @ARGV;

my @fasta;
my $outdir = 'ntBiases';
my $minsize = 1;
my $winsize = 10000;
my $step = 5000;
my $reference;
my $ncheck;
my $gap = 0;
my $unit = 'Mb';
my $labels = 'mixed';
my $label_size = 36;
my $label_font = 'bold';
my $circos_plot;
my $custom_file;
my $custom_cc;
my $list_preset;
my $tsv;
my $noticks;
my $max_ticks = 5000;
my $max_ideograms = 200;
my $max_links = 75000;
my $max_points_per_track = 75000;
my $zdepth = 50;
my $clusters;
my $no_biases;
GetOptions(
	'f|fasta=s@{1,}' => \@fasta,
	'o|outdir=s' => \$outdir,
	'minsize=i' => \$minsize,
	'w|winsize=i' => \$winsize,
	's|step=i' => \$step,
	'n|ncheck' => \$ncheck,
	'c|circos' => \$circos_plot,
	'u|unit=s' => \$unit,
	'l|labels=s' => \$labels,
	'label_size=s' => \$label_size,
	'label_font=s' => \$label_font,
	'r|reference=s' => \$reference,
	'g|gap=i' => \$gap,
	'custom_file=s' => \$custom_file,
	'custom_preset=s' => \$custom_cc,
	'list_preset'	=> \$list_preset,
	't|tsv' => \$tsv,
	'noticks' => \$noticks,
	'max_ticks=i' => \$max_ticks,
	'max_ideograms=i' => \$max_ideograms,
	'max_links=i' => \$max_links,
	'max_points_per_track=i' => \$max_points_per_track,
	'zdepth=i' => \$zdepth,
	'clusters' => \$clusters,
	'no_biases' => \$no_biases
);

### List presets and stop
if ($list_preset){

	print "\n".'Available color presets:'."\n";

	my %available_colors = cc_colors();
	foreach my $color_key (sort (keys %available_colors)){
		my $color_number = scalar (keys %{$available_colors{$color_key}});
		print $color_key."\t".$color_number." colors\n";
	}

	print "\n";
	exit;

}

### Check if output directory / subdirs can be created
$outdir =~ s/\/$//;
my $catdir = $outdir.'/concatenated';
my $pairwisedir = $outdir.'/pairwise';
my $singledir = $outdir.'/single';

for my $dir ($outdir,$catdir,$pairwisedir,$singledir){
	unless (-d $dir) {
		make_path($dir,{mode => 0755}) or die "Can't create $dir: $!\n";
	}
}

no warnings 'once';
my @cathandles = (*CGC, *CAT, *CGT, *CAC, *CGA, *CCT);
foreach my $cfh (@cathandles){
	my ($lfh) = $cfh =~ /C(\w+)$/;
	my $cat_bias = $catdir.'/concatenated.'.$lfh;
	open $cfh, '>', $cat_bias or die "Can't create $cat_bias: $!\n";
	print $cfh '#chr START END GC_PERCENTAGE'."\n";
}

### Iterating through FASTA file(s); creating database of sequences (could be multifasta)

my %sequences;
my %percent;
my $seqname;
my $fileprefix;

while (my $fasta = shift@fasta){

	my ($basename) = fileparse($fasta);
	($fileprefix) = $basename =~ /(\S+)\.\w+$/;

	open FASTA, "<", $fasta or die "Can't open $fasta: $!\n";

	## Database of sequences
	while (my $line = <FASTA>){
		chomp $line;
		if ($line =~ /^>(\S+)/){
			$seqname = $1;
		}
		else {
			$sequences{$fileprefix}{$seqname} .= $line;
		}
	}

	## Delete sequence if smaller than minsize
	foreach my $seq (keys %{$sequences{$fileprefix}}){
		my $seq_len = length($sequences{$fileprefix}{$seq});
		if ($seq_len < $minsize){
			delete $sequences{$fileprefix};
		}
	}

	my $fasta_dir = $singledir.'/'.$fileprefix;
	unless (-d $fasta_dir){
		mkdir ($fasta_dir,0755) or die "Can't create $fasta_dir: $!\n";
	}

	### Circos
	no warnings 'once'; ## Filehandles are used in a loop...
	my @filehandles = (*GC, *AT, *GT, *AC, *GA, *CT);

	if ($ncheck){
		@filehandles = (*GC, *AT, *GT, *AC, *GA, *CT, *NN);
	}

	my $circos_GC;
	my $circos_AT;
	my $circos_GT;
	my $circos_AC;
	my $circos_GA;
	my $circos_CT;
	my $circos_NN;

	## Nucleotide biases
	for my $fh (@filehandles){
		my ($lfh) = $fh =~ /(\w+)$/; ## grabbing test from filehandle
		my $subdir = $singledir.'/'.$fileprefix;
		my $filename = $subdir.'/'.$fileprefix.'.'.$lfh;
		open $fh, ">", $filename or die "Can't create $filename: $!\n";
		print $fh '#chr START END GC_PERCENTAGE'."\n";
	}

	### Iterating through each sequence in the FASTA file
	foreach my $sequence (sort (keys %{$sequences{$fileprefix}})){

		my $outfile = $fasta_dir.'/'.$fileprefix.'.'.$sequence.'.tsv';

		if ($tsv){
			open BIAS, ">", $outfile or die "Can't create $outfile: $!\n";
			print BIAS "# Location\t% GC\t% AT\t% AG\t% CT\t% GT\t% AC";
			if ($ncheck){
				print BIAS "\t% NN";
			}
			print BIAS "\n";
		}

		### Sliding windows
		my $seq = $sequences{$fileprefix}{$sequence};
		my $csize = length $seq;
		my $x;

		for ($x = 0; $x <= ($csize - $winsize); $x += $step){

			my $subseq = substr($seq, $x, $winsize);
			my $end = $x + $winsize - 1;
			
			%percent = ();
			biases($subseq,$x,$winsize);

			foreach my $fh (@filehandles){
				my ($lfh) = $fh =~ /(\w+)$/;
				print $fh "$sequence $x $end $percent{$lfh}\n";
			}
			foreach my $cfh (@cathandles){
				my ($lfh) = $cfh =~ /C(\w+)$/;
				print $cfh "$sequence $x $end $percent{$lfh}\n";
			}

		}

		### Working on leftover string < $winsize
		my $modulo = $csize % $winsize;
		my $subseqleft = substr($seq, -$modulo, $modulo);
		my $leftover_size = length $subseqleft;

		%percent = ();
		biases($subseqleft,$x,$leftover_size);

		my $end = $csize - 1;
		foreach my $fh (@filehandles){
			my ($lfh) = $fh =~ /(\w+)$/;
			print $fh "$sequence $x $end $percent{$lfh}\n";
		}
		foreach my $cfh (@cathandles){
			my ($lfh) = $cfh =~ /C(\w+)$/;
			print $cfh "$sequence $x $end $percent{$lfh}\n";
		}

		close BIAS;

	}

	close FASTA;

}

####################################################################
### Circos
####################################################################

############## Creating karyotype files

# Normal = karyotypes in the same order as encountered
# Inverted = karyotypes in reversed order (except in the contatenated file):
# In the inverted concatenated file, all but the reference karyotypes
# are in reversed order. This can help visualize synteny when comparing
# genomes to the reference => sometimes the reversed order is more
# informative than the normal one

my $cat_kar = $catdir.'/concatenated.normal.karyotype';
my $cat_kar_inv = $catdir.'/concatenated.inverted.karyotype';

open CONCAT, ">", $cat_kar or die "Can't create $cat_kar: $!\n";
open CONCATINV, ">", $cat_kar_inv or die "Can't create $cat_kar_inv: $!\n";

print CONCAT '#chr - ID LABEL START END COLOR'."\n";
print CONCATINV '#chr - ID LABEL START END COLOR'."\n";

## Using first sequence as reference if none is provided
	if (!$reference){
		my @genomes = sort (keys %sequences);
		$reference = $genomes[0];
	}

foreach my $fileprefix (keys (%sequences)){

	## Creating a "karyotype" file for Circos
	my $circos_dir = $singledir.'/'.$fileprefix;
	my $circos_kar = $circos_dir.'/'.$fileprefix.'.normal.karyotype';
	my $circos_kar_inv = $circos_dir.'/'.$fileprefix.'.inverted.karyotype';

	open KAR, ">", $circos_kar or die "Can't create $circos_kar: $!\n";
	open KARINV, ">", $circos_kar_inv or die "Can't create $circos_kar_inv: $!\n";

	print KAR '#chr - ID LABEL START END COLOR'."\n";
	print KARINV '#chr - ID LABEL START END COLOR'."\n";

	my @seqs = sort (keys %{$sequences{$fileprefix}});
	my $num_of_seq = scalar (@seqs);

	my $id = 0;
	my $label;

	foreach my $sequence (@seqs){

		my $seq = $sequences{$fileprefix}{$sequence};
		my $csize = length $seq;

		my $terminus = $csize - 1;
		$id++;

		if ($labels eq 'names'){
			$label = $sequence;
		}
		else {
			$label = $id;
		}

		if ($labels eq 'mixed'){
			if ($fileprefix eq $reference){
				$label = uc(int2roman($label));
			}
		}
		elsif ($labels eq 'roman'){
			$label = uc(int2roman($label));
		}

		print KAR "chr - $sequence $label 0 $terminus black\n";
		print CONCAT "chr - $sequence $label 0 $terminus black\n";

		if ($fileprefix eq $reference){
			print CONCATINV "chr - $sequence $label 0 $terminus black\n";
		}

	}

	my @rev_seqs = reverse(@seqs);
	my $id_rev = $num_of_seq;

	foreach my $sequence (@rev_seqs){

		my $seq = $sequences{$fileprefix}{$sequence};
		my $csize = length $seq;

		my $terminus = $csize - 1;

		if ($labels eq 'names'){
			$label = $sequence;
		}
		else {
			$label = $id_rev;
		}

		if ($labels eq 'mixed'){
			if ($fileprefix eq $reference){
				$label = uc(int2roman($label));
			}
		}
		elsif ($labels eq 'roman'){
			$label = uc(int2roman($label));
		}

		print KARINV "chr - $sequence $label 0 $terminus black\n";

		if ($fileprefix ne $reference){
			print CONCATINV "chr - $sequence $label 0 $terminus black\n";
		}
		
		$id_rev--;

	}

	close KAR;
	close KARINV;

}

close CONCAT;
close CONCATINV;

### Pairwise
foreach my $query (keys (%sequences)){

	foreach my $subject (keys (%sequences)){

		my $pairwise_ref = $query; ## Setting the query as the reference in pairwise mode

		if ($query eq $subject){
			next;
		}
		else {

			## Creating a "karyotype" file for Circos
			my $prefix = $query.'_vs_'.$subject;
			my $subdir = $pairwisedir.'/'.$prefix;
			unless (-d $subdir){
				mkdir($subdir,0755) or die "Can't create $subdir: $!\n";
			}
			my $circos_kar = $subdir.'/'.$prefix.'.normal.karyotype';
			my $circos_kar_inv = $subdir.'/'.$prefix.'.inverted.karyotype';

			open PAIR, ">", $circos_kar or die "Can't create $circos_kar: $!\n";
			open PAIRINV, ">", $circos_kar_inv or die "Can't create $circos_kar_inv: $!\n";

			print PAIR '#chr - ID LABEL START END COLOR'."\n";
			print PAIRINV '#chr - ID LABEL START END COLOR'."\n";

			for my $key ($query, $subject){

				my @seqs = sort (keys %{$sequences{$key}});
				my $num_of_seq = scalar (@seqs);
				my $id = 0;
				my $label;

				foreach my $sequence (@seqs){

					my $seq = $sequences{$key}{$sequence};
					my $csize = length $seq;

					my $terminus = $csize - 1;
					$id++;

					if ($labels eq 'names'){
						$label = $sequence;
					}
					else {
						$label = $id;
					}

					if ($labels eq 'mixed'){
						if ($key eq $pairwise_ref){
							$label = uc(int2roman($label));
						}
					}
					elsif ($labels eq 'roman'){
						$label = uc(int2roman($label));
					}

					print PAIR "chr - $sequence $label 0 $terminus black\n";

					if ($key eq $pairwise_ref){
						print PAIRINV "chr - $sequence $label 0 $terminus black\n";
					}

				}

				my @rev_seqs = reverse(@seqs);
				my $id_rev = $num_of_seq;

				foreach my $sequence (@rev_seqs){

					my $seq = $sequences{$key}{$sequence};
					my $csize = length $seq;

					my $terminus = $csize - 1;

					if ($labels eq 'names'){
						$label = $sequence;
					}
					else {
						$label = $id_rev;
					}

					if ($labels eq 'mixed'){
						if ($key eq $pairwise_ref){
							$label = uc(int2roman($label));
						}
					}
					elsif ($labels eq 'roman'){
						$label = uc(int2roman($label));
					}

					if ($key ne $pairwise_ref){
						print PAIRINV "chr - $sequence $label 0 $terminus black\n";
					}
					
					$id_rev--;

				}

			}

			close PAIR;
			close PAIRINV;

		}

	}

}


############## Creating a default ideogram configuration for Circos
my $ideogram = $outdir.'/'.'ideogram.conf';
open my $id, '>', $ideogram or die $!;

my $ideogram_data = <<"IDEO";
<ideogram>

<spacing>
default = 0.01r
</spacing>

radius    = 0.9r
thickness = 5p
fill      = yes
show_label       = yes
label_font       = $label_font
label_radius     = dims(ideogram,radius) + 0.07r
label_size       = $label_size
label_parallel   = yes

</ideogram>
IDEO
print $id $ideogram_data."\n";
close $id;

############## Creating a default ticks configuration for Circos
my $ticks = $outdir.'/'.'ticks.conf';
open my $tk, '>', $ticks or die $!;

my $multiplier = '1e-6';
if ($unit =~ /kb/i){
	$multiplier = '1e-3';
}

my $ticks_data = <<"TICKS";
show_ticks          = yes
show_tick_labels    = yes

<ticks>
radius           = 1r
orientation	= out
tick_separation  = 3p
label_separation = 1p
multiplier       = $multiplier
color            = black
thickness        = 3p

<tick>
size             = 20p
spacing        = 10u
show_label     = yes
label_size     = 15p
label_offset   = 5p
format         = \%d $unit
grid           = yes
grid_color     = dgrey
grid_thickness = 1p
grid_start     = 1r
grid_end       = 1.1r
</tick>

<tick>
size             = 10p
spacing        = 5u
show_label     = yes
label_size     = 10p
label_offset   = 5p
format         = \%.1f
thickness      = 1.5p
color          = vdgrey
</tick>

<tick>
size             = 5p
spacing        = 1u
show_label     = no
thickness      = 1p
color          = dgrey
</tick>

</ticks>
TICKS
print $tk $ticks_data."\n";
close $tk;

#################### Circos colors

my @color_set;

# 6 colors per set
my @yellows = ('vlyellow','lyellow','yellow','dyellow','vdyellow','vvdyellow');

# 7 colors per set
my @reds = ('vvlred','vlred','lred','red','dred','vdred','vvdred');
my @oranges = ('vvlorange','vlorange','lorange','orange','dorange','vdorange','vvdorange');
my @greens = ('vvlgreen','vlgreen','lgreen','green','dgreen','vdgreen','vvdgreen');
my @blues = ('vvlblue','vlblue','lblue','blue','dblue','vdblue','vvdblue');
my @purples = ('vvlpurple','vlpurple','lpurple','purple','dpurple','vdpurple','vvdpurple');

# 8 colors per set; did not include vvvlgrey => too light can't see it on a white background
my @greys = ('vvlgrey','vlgrey','lgrey','grey','dgrey','vdgrey','vvdgrey','vvvdgrey');

my @rainbow = (@reds,@oranges,@yellows,@greens,@blues,@purples); ## 41 colors total
my @bowgrey = (@rainbow,@greys); ## 49 colors total

# Custom color set(s)
my %custom_colors = cc_colors();
my @custom_set;

if ($custom_file){

	open CC, '<', $custom_file or die "Can't read $custom_file: $!\n";
	my ($basename) = fileparse($custom_file);
	$custom_cc = $basename;

	while (my $line = <CC>){
		chomp $line;
		unless ($line =~ /^#/){
			my ($color,$rgb) = split("\t", $line);
			$custom_colors{$basename}{$color} = $rgb;
		}

	}

}

if ($custom_cc){
	@custom_set = sort (keys %{$custom_colors{$custom_cc}});
} 

############## Creating configuration files for Circos

### Nucleotides biases calculated by nucleotide_biases.pl
### Using same color scheme as in Mascarenhas dos Santos et al. (2023)
### https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-023-09331-3

my @biases = ('GC','AT','GT','AC','GA','CT');
if ($ncheck){
	@biases = ('GC','AT','GT','AC','GA','CT','NN');
}
my %colors = (
	'GC' => 'vdred',
	'AT' => 'vdgrey',
	'GT' => 'vdblue',
	'AC' => 'vdgreen',
	'GA' => 'vdpurple',
	'CT' => 'vdyellow',
	'NN' => 'black'
);

### Creating files for both normal and inverted karyotypes (single and concatenated)
for my $genome ((keys %sequences), 'concatenated'){

	my $subfile;
	my $subdir = $singledir.'/'.$genome;

	if ($genome eq 'concatenated'){
		$subfile = $outdir.'/concatenated/'.$genome;
	}
	else{
		$subfile = $subdir.'/'.$genome;
	}

	my $pairwise_ref = $reference;

	## Creating conf files for Circos
	for my $orientation ('normal', 'inverted'){
		if ($genome eq 'concatenated'){
			print_circos_conf($subfile,$orientation,$pairwise_ref,'cat','links');
		}
		else{
			print_circos_conf($subfile,$orientation,$pairwise_ref,'nocat','nolink');
		}
	}

}

### Creating pairwise files for both normal and inverted karyotypes
for my $genome (keys %sequences){

	for my $subject (keys %sequences){

		if ($genome eq $subject){
			next;
		}

		## Creating conf files for Circos
		my $subdir = $pairwisedir.'/'.$genome.'_vs_'.$subject;
		my $subfile = $subdir.'/'.$genome.'_vs_'.$subject;
		my $pairwise_ref = $genome; ## Setting the query as the reference in pairwise mode

		for my $orientation ('normal', 'inverted'){
			print_circos_conf($subfile,$orientation,$pairwise_ref,'cat','links');
		}

	}

}

####################################################################
### Subroutine(s)
####################################################################

sub print_circos_conf {

	my $subfile = $_[0];
	my $orientation = $_[1];
	my $pairwise_ref = $_[2];
	my $cat_status = $_[3];
	my $link_status = $_[3];

	my $config = $subfile.'.'.$orientation.'.conf';
	my $karyotype = $subfile.'.'.$orientation.'.karyotype';

	my $bias_file = $subfile;
	if ($cat_status){
		$bias_file = $catdir.'/concatenated';
	}

	open my $cg, '>', $config or die "Can't create $config: $!\n";

		print $cg "<<include $ideogram>>"."\n";
		if ($noticks){
			print $cg '#';
		}
		print $cg "<<include $ticks>>"."\n\n";
		print $cg "karyotype = $karyotype"."\n";

		my $unit_width = 100000;
		if ($unit =~ /kb/i){
			$unit_width = 10000;
		}
		print $cg 'chromosomes_units='.$unit_width."\n\n";



		## Biases plots
		my $r_start = 0.99;
		my $r_end = 0.95;

		if ($no_biases){
			goto LINKS;
		}

		print $cg '<plots>'."\n";
		print $cg '################# BIASES'."\n\n";

		my $modulo_counter = 0;

		for my $bias (@biases) {

			$modulo_counter++;

			print $cg '## '.$bias."\n";
			print $cg '<plot>'."\n";
			print $cg "file    = $bias_file.$bias"."\n";
			print $cg 'type      = line'."\n";
			print $cg 'thickness = 2'."\n";
			print $cg 'max_gap = 1u'."\n";
			print $cg "color   = ".$colors{$bias}."\n";
			print $cg 'min     = 20'."\n";
			print $cg 'max     = 80'."\n";
			print $cg "r1      = ${r_start}r"."\n";
			print $cg "r0      = ${r_end}r"."\n\n";
			axes(\*$cg);
			print $cg '</plot>'."\n\n";

			if (($modulo_counter % 2) == 0){
				$r_start -= 0.05;
				$r_end -= 0.05;
			}

		}

		print $cg '</plots>'."\n\n";

		LINKS:
		unless ($link_status eq 'nolink'){

			my $link_start = $r_end + 0.04;
			my $link_color = 'grey_a5';

			print $cg '########### Links'."\n\n";
			print $cg '<links>'."\n\n";

			print $cg '<link>'."\n";
			print $cg "file          = $subfile.gap_$gap.links"."\n";
			print $cg "radius        = ${link_start}r"."\n";
			print $cg 'bezier_radius = 0.1r'."\n";
			print $cg 'ribbon = yes'."\n";
			print $cg 'color         = '.$link_color."\n";
			print $cg 'thickness     = 1'."\n";

			unless ($clusters){

				## Rules
				print $cg '<rules>'."\n\n";

				## Use first entry if pairwise reference does not exits
				unless (exists $sequences{$pairwise_ref}){
					my @tmp = sort (keys %sequences);
					$pairwise_ref = $tmp[0];
				}

				## Counting for required colors
				my $ref_sequence_count = scalar (keys %{$sequences{$pairwise_ref}});

				if ($custom_cc){
					@color_set = @custom_set;
				}
				elsif ($ref_sequence_count <= scalar(@rainbow)){
					@color_set = @rainbow;
				}
				else {
					@color_set = @bowgrey;
				}

				## Creating an increment so that it will use the full range of colors
				## not just the start of the list
				my $increment = scalar(@color_set)/$ref_sequence_count;
				my $rounded_increment = int(sprintf("%.0f", $increment));

				## Making sure that the increment is >= 1
				if ($rounded_increment <= 1){
					$rounded_increment = 1;
				}

				my $color_start = 0;
				
				foreach my $refseq (sort (keys %{$sequences{$pairwise_ref}})){
					foreach my $queseq (sort (keys %sequences)){
						if (($queseq eq 'concatenated') or ($queseq eq $pairwise_ref)){
							next;
						}
						else{
							foreach my $queseq_cg (sort (keys %{$sequences{$queseq}})){
								print $cg '<rule>'."\n";
								print $cg "condition  = between($refseq,$queseq_cg)"."\n";
								print $cg "color      = ".$color_set[$color_start]."\n";
								print $cg 'flow       = continue'."\n";
								print $cg 'z          = '.$zdepth."\n";
								print $cg '</rule>'."\n\n";
							}
						}
					}
					$color_start += $rounded_increment;
					## Check if no more colors left, if so restart from 1st color in color set
					if ($color_start >= scalar(@color_set)){
						$color_start = 0;
					}
				}
				print $cg '</rules>'."\n\n";

			}

			print $cg '</link>'."\n";
			print $cg '</links>'."\n\n";

		}

		## image.conf
		print $cg '### Image conf:'."\n";
		print $cg '<image>'."\n";
		image_conf(\*$cg);
		print $cg '</image>'."\n\n";

		## colors.conf
		print $cg '<<include etc/colors_fonts_patterns.conf>>'."\n";

		## housekeeping
		print $cg "\n".'### Housekeeping:'."\n";
		housekeeping(\*$cg);

		## add custom colors (if requested)
		if ($custom_cc){
			print $cg '<colors>'."\n";
			foreach my $color (@custom_set){
				print $cg $color.' = '.$custom_colors{$custom_cc}{$color}."\n";
			}
			print $cg '</colors>'."\n";
		}

		close $cg;

}

sub biases {

		my $curseq = $_[0];
		my $pos = $_[1];
		my $divider = $_[2];

		if ($divider == 0){
			next;
		}

		my $gc = $curseq =~ tr/GgCc//;
		my $at = $curseq =~ tr/AaTt//;
		my $ga = $curseq =~ tr/GgAa//;
		my $ct = $curseq =~ tr/CcTt//;
		my $gt = $curseq =~ tr/GgTt//;
		my $ac = $curseq =~ tr/AaCc//;
		my $nn = $curseq =~ tr/Nn//;

		$percent{'GC'} = $gc = ($gc/$divider) * 100;
		$percent{'AT'} = $at = ($at/$divider) * 100;
		$percent{'GA'} = $ga = ($ga/$divider) * 100;
		$percent{'CT'} = $ct = ($ct/$divider) * 100;
		$percent{'GT'} = $gt = ($gt/$divider) * 100;
		$percent{'AC'} = $ac = ($ac/$divider) * 100;
		$percent{'NN'} = $nn = ($nn/$divider) * 100;

		if ($tsv){
			print BIAS "$pos\t$gc\t$at\t$ga\t$ct\t$gt\t$ac";
			if ($ncheck){
				print BIAS "\t$nn";
			}
			print BIAS "\n";
		}

}

sub axes {
	my $fh = $_[0];
	print $fh '<axes>'."\n";
	print $fh "\t".'<axis>'."\n";
	print $fh "\t".'color     = lgrey_a2'."\n";
	print $fh "\t".'thickness = 1'."\n";
	print $fh "\t".'spacing   = 0.025r'."\n";
	print $fh "\t".'</axis>'."\n";
	print $fh '</axes>'."\n";
}

sub image_conf {
	my $fh = $_[0];
	print $fh 'background = '.'white'."\n";
	print $fh 'dir   = . '."\n";
	print $fh 'file  = '.'circos.png'."\n";
	print $fh 'png   = '.'yes'."\n";
	print $fh 'svg   = '.'yes'."\n";
	print $fh 'radius         = '.'1500p'."\n";
	print $fh 'angle_offset      = '.'-90'."\n";
	print $fh 'auto_alpha_colors = '.'yes'."\n";
	print $fh 'auto_alpha_steps  = '.'5'."\n";
}

sub cc_colors {
	%custom_colors = (
		'chloropicon' => { # from Lemieux et al. (2019) https://www.nature.com/articles/s41467-019-12014-x
			'c01' => '202,75,75',
			'c02' => '239,60,104',
			'c03' => '241,102,140',
			'c04' => '245,152,162',
			'c05' => '245,126,47',
			'c06' => '250,166,55',
			'c07' => '255,197,61',
			'c08' => '255,228,67',
			'c09' => '210,213,76',
			'c10' => '147,195,84',
			'c11' => '18,178,89',
			'c12' => '0,179,127',
			'c13' => '0,180,161',
			'c14' => '0,182,204',
			'c15' => '0,183,241',
			'c16' => '0,157,218',
			'c17' => '64,131,196',
			'c18' => '94,104,176',
			'c19' => '108,82,162',
			'c20' => '122,42,144'
		},
		'encephalitozoon' => { # from Pombert et al. (2012) https://pubmed.ncbi.nlm.nih.gov/22802648/
			'ei-01' => '0,162,120',
			'ei-02' => '0,165,180',
			'ei-03' => '91,202,244',
			'ei-04' => '139,162,211',
			'ei-05' => '121,97,169',
			'ei-06' => '162,25,141',
			'ei-07' => '235,0,139',
			'ei-08' => '240,102,129',
			'ei-09' => '241,101,80',
			'ei-10' => '245,138,32',
			'ei-11' => '164,115,11'
		},
		'blues' => { # A simple blue gradient (13 hues)
			'bl-01' => '245,255,255',
			'bl-02' => '190,233,244',
			'bl-03' => '180,212,233',
			'bl-04' => '170,191,222',
			'bl-05' => '160,171,211',
			'bl-06' => '149,150,200',
			'bl-07' => '138,130,190',
			'bl-08' => '127,111,179',
			'bl-09' => '115,91,168',
			'bl-10' => '103,71,157',
			'bl-11' => '91,52,146',
			'bl-12' => '77,30,136',
			'bl-13' => '63,0,125'
		},
		### add custom colors as desired
		# 'custom_name_2' => {
		# 	'key_1' => 'rgb_value_1',
		# 	'key_2' => 'rgb_value_2',
		# 	'key_3' => 'rgb_value_3',
		# },
		# 'custom_name_2' => {
		# 	'key_1' => 'rgb_value_1',
		# 	'key_2' => 'rgb_value_2',
		# 	'key_3' => 'rgb_value_3',
		# },
	)
}

sub housekeeping {
	
	my $fh = $_[0]; 

	print $fh 'anglestep       = '.'0.5'."\n";
	print $fh 'minslicestep    = '.'10'."\n";
	print $fh 'beziersamples   = '.'40'."\n";
	print $fh 'debug           = '.'no'."\n";
	print $fh 'warnings        = '.'no'."\n";
	print $fh 'imagemap        = '.'no'."\n";
	print $fh 'paranoid        = '.'yes'."\n";
	print $fh 'units_ok        = '.'bupr'."\n";
	print $fh 'units_nounit    = '.'n'."\n";
	print $fh 'file_delim = '.'\s'."\n";
	print $fh 'file_delim_collapse = '.'yes'."\n";
	print $fh 'list_record_delim = '.'\s*[;,]\s*'."\n";
	print $fh 'list_field_delim  = '.'\s*[:=]\s*'."\n";
	print $fh 'options_record_delim = '.'[,;]'."\n";
	print $fh 'options_field_delim  = '.'='."\n";
	print $fh 'skip_missing_expression_vars = '.'no'."\n";
	print $fh 'legacy_underline_expression_syntax = '.'no'."\n";
	print $fh 'svg_font_scale = '.'1.3'."\n";
	print $fh 'sup_baseline_shift = '.'40'."\n";
	print $fh 'sub_baseline_shift = '.'-40'."\n";
	print $fh 'sup_fontsize = '.'90'."\n";
	print $fh 'sub_fontsize = '.'90'."\n";
	print $fh 'default_font   = '.'default'."\n";
	print $fh 'default_font_name  = '.'Arial'."\n";
	print $fh 'default_font_color = '.'black'."\n";
	print $fh 'default_color  = '.'black'."\n";
	print $fh '<guides>'."\n";
	print $fh 'thickness      = '.'1'."\n";
	print $fh 'size           = '.'5'."\n";
	print $fh 'type           = '.'outline'."\n";
	print $fh '<object>'."\n";
	print $fh 'all            = '.'no'."\n";
	print $fh 'ideogram       = '.'no'."\n";
	print $fh 'ideogram_label = '.'no'."\n";
	print $fh '</object>'."\n";
	print $fh '<color>'."\n";
	print $fh 'default = '.'lblue'."\n";
	print $fh 'text    = '.'red'."\n";
	print $fh '</color>'."\n";
	print $fh '</guides>'."\n";
	print $fh 'debug_group = '.'summary,output'."\n";
	print $fh 'debug_auto_timer_report = '.'30'."\n";
	print $fh 'debug_word_separator = '.'" "'."\n";
	print $fh 'debug_undef_text     = '.'_undef_'."\n";
	print $fh 'debug_empty_text     = '.'_emptylist_'."\n";
	print $fh 'debug_validate       = '.'yes'."\n";
	print $fh 'debug_output_tidy    = '.'no'."\n";
	print $fh 'text_pixel_subsampling = '.'1'."\n";
	print $fh 'text_snuggle_method    = '.'array'."\n";
	print $fh 'restrict_parameter_names = '.'no'."\n";
	print $fh 'case_sensitive_parameter_names = '.'no'."\n";
	print $fh 'calculate_track_statistics = '.'yes'."\n";
	print $fh 'color_cache_static = '.'yes'."\n";
	print $fh 'color_cache_file   = '.'circos.colorlist'."\n";
	print $fh 'color_lists_use    = '.'yes'."\n";
	print $fh 'memoize = '.'yes'."\n";
	print $fh 'quit_on_dump = '.'yes'."\n";
	print $fh 'offsets = '.'0,0'."\n";
	print $fh 'max_ticks            = '.$max_ticks."\n";
	print $fh 'max_ideograms        = '.$max_ideograms."\n";
	print $fh 'max_links            = '.$max_links."\n";
	print $fh 'max_points_per_track = '.$max_points_per_track."\n";
	print $fh 'undefined_ideogram = '.'skip'."\n";
	print $fh 'relative_scale_iterations = '.'10'."\n";
	print $fh 'relative_scale_spacing    = '.'mode'."\n";
	print $fh 'data_out_of_range = '.'trim,warn'."\n";
	print $fh 'track_defaults = '.'etc/tracks'."\n";
	print $fh 'round_brush_use           = '.'yes'."\n";
	print $fh 'round_brush_min_thickness = '.'5'."\n";
	print $fh 'anti_aliasing = '.'yes'."\n";
	print $fh 'housekeeping = '.'yes'."\n";
	print $fh 'auto_eval = '.'no'."\n";

}