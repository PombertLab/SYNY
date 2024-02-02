#!/usr/bin/perl
## Pombert Lab, 2022
my $name = 'nucleotide_biases.pl';
my $version = '0.4f';
my $updated = '2024-02-02';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;
use File::Path qw(make_path);
use Math::Round;

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
		  -reference CCMP1205 \\
		  -circos \\
		  -gap 0 \\
		  -custom

-f (--fasta)	Fasta file(s) to process
-o (--outdir)	Output directory [Default: ntBiases]
-w (--winsize)	Sliding window size [Default: 10000]
-s (--step)		Sliding window step [Default: 5000]
-r (--reference)	Genome reference for Circos plotting
-n (--ncheck)	Check for ambiguous/masked (Nn) nucleotides
-c (--circos)	Run Circos to plot images
-g (--gap)		Default gap links file for Circos plotting [Default: 0]
-u (--unit)		Size unit (Kb or Mb) [Default: Mb]
-custom		Use custom colors [c01 to c20]
-t (--tsv)		Output tab-delimited files (e.g. for excel plotting)
OPTIONS
die "\n$usage\n" unless @ARGV;

my @fasta;
my $outdir = 'ntBiases';
my $winsize = 1000;
my $step = 500;
my $reference;
my $ncheck;
my $gap = 0;
my $unit = 'Mb';
my $circos_plot;
my $custom_cc;
my $tsv;
GetOptions(
	'f|fasta=s@{1,}' => \@fasta,
	'o|outdir=s' => \$outdir,
	'w|winsize=i' => \$winsize,
	's|step=i' => \$step,
	'n|ncheck' => \$ncheck,
	'c|circos' => \$circos_plot,
	'u|unit=s' => \$unit,
	'r|reference=s' => \$reference,
	'g|gap=i' => \$gap,
	'custom' => \$custom_cc,
	't|tsv' => \$tsv
);

### Check if output directory / subdirs can be created
$outdir =~ s/\/$//;
unless (-d $outdir) {
	make_path($outdir,{mode => 0755}) or die "Can't create $outdir: $!\n";
}
my $catdir = $outdir.'/CONCATENATED';
unless (-d $catdir) {
	make_path($catdir,{mode => 0755}) or die "Can't create $catdir: $!\n";
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

	while (my $line = <FASTA>){
		chomp $line;
		if ($line =~ /^>(\S+)/){
			$seqname = $1;
		}
		else {
			$sequences{$fileprefix}{$seqname} .= $line;
		}
	}

	my $fasta_dir = "$outdir/$fileprefix";
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
	
	my $circos_dir = $outdir.'/'.$fileprefix;
	unless (-d $circos_dir) {
			make_path($circos_dir,{mode => 0755}) or die "Can't create $circos_dir: $!\n";
	}

	## Nucleotide biases
	for my $fh (@filehandles){
		my ($lfh) = $fh =~ /(\w+)$/; ## grabbing test from filehandle
		my $filename = $circos_dir.'/'.$fileprefix.'.'.$lfh;
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

############## Creating genotype/karyotype files

# Normal = genotypes in the same order as encountered
# Inverted = genotypes in reversed order (except in the contatenated file):
# In the inverted concatenated file, all but the reference genotypes
# are in reversed order. This can help visualize synteny when comparing
# genomes to the reference => sometimes the reversed order is more
# informative than the normal one

my $cat_kar = $catdir.'/concatenated.genotype';
my $cat_kar_inv = $catdir.'/concatenated.inverted.genotype';

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
	my $circos_dir = $outdir.'/'.$fileprefix;
	my $circos_kar = $circos_dir.'/'.$fileprefix.'.genotype';
	my $circos_kar_inv = $circos_dir.'/'.$fileprefix.'.inverted.genotype';

	open KAR, ">", $circos_kar or die "Can't create $circos_kar: $!\n";
	open KARINV, ">", $circos_kar_inv or die "Can't create $circos_kar_inv: $!\n";

	print KAR '#chr - ID LABEL START END COLOR'."\n";
	print KARINV '#chr - ID LABEL START END COLOR'."\n";

	my @seqs = sort (keys %{$sequences{$fileprefix}});
	my $num_of_seq = scalar (@seqs);

	my $id = 0;

	foreach my $sequence (@seqs){

		my $seq = $sequences{$fileprefix}{$sequence};
		my $csize = length $seq;

		my $terminus = $csize - 1;
		$id++;
		print KAR "chr - $sequence $id 0 $terminus black\n";
		print CONCAT "chr - $sequence $id 0 $terminus black\n";

		if ($fileprefix eq $reference){
			print CONCATINV "chr - $sequence $id 0 $terminus black\n";
		}

	}

	
	my @rev_seqs = reverse(@seqs);
	my $id_rev = $num_of_seq;

	foreach my $sequence (@rev_seqs){

		my $seq = $sequences{$fileprefix}{$sequence};
		my $csize = length $seq;

		my $terminus = $csize - 1;
		print KARINV "chr - $sequence $id_rev 0 $terminus black\n";

		if ($fileprefix ne $reference){
			print CONCATINV "chr - $sequence $id_rev 0 $terminus black\n";
		}
		
		$id_rev--;

	}

	close KAR;
	close KARINV;

}

close CONCAT;
close CONCATINV;

############## Creating a default ideogram configuration for Circos
my $ideogram = $outdir.'/'.'ideogram.conf';
open my $id, '>', $ideogram or die $!;

my $ideogram_data = <<'IDEO';
<ideogram>

<spacing>
default = 0.01r
</spacing>

radius    = 0.9r
thickness = 5p
fill      = yes
show_label       = yes
label_font       = bold 
label_radius     = dims(ideogram,radius) + 0.07r
label_size       = 36
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

# 20 custom color set 
my @custom_set = (
	'c01','c02','c03','c04','c05','c06','c07','c08','c09','c10',
	'c11','c12','c13','c14','c15','c16','c17','c18','c19','c20'
);

custom_colors();


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

### Creating files for both normal and inverted genotypes
for my $genome ((keys %sequences), 'concatenated'){

	my $subdir = $outdir.'/'.$genome;
	my $subfile = $subdir.'/'.$genome;
	my $karyotype = $subfile.'.genotype';
	my $config = $subfile.'.conf';
	
	## Creating conf files for Circos
	for my $orientation ('normal', 'inverted'){

		if ($orientation eq 'inverted'){
			$karyotype = $subfile.'.inverted.genotype';
			$config = $subfile.'.inverted.conf';
		}

		open my $cg, '>', $config or die $!;

		print $cg "<<include $ideogram>>"."\n";
		print $cg "<<include $ticks>>"."\n\n";
		print $cg "karyotype = $karyotype"."\n";

		my $unit_width = 100000;
		if ($unit =~ /kb/i){
			$unit_width = 10000;
		}
		print $cg 'chromosomes_units='.$unit_width."\n\n";

		## Biases plots
		print $cg '<plots>'."\n";
		print $cg '################# BIASES'."\n\n";

		my $r_start = 0.99;
		my $r_end = 0.95;
		my $modulo_counter = 0;

		for my $bias (@biases) {

			$modulo_counter++;

			print $cg '## '.$bias."\n";
			print $cg '<plot>'."\n";
			print $cg "file    = $subfile.$bias"."\n";
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

		## Links; only for concatenated genomes
		if ($genome eq 'concatenated'){
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

			## Rules
			print $cg '<rules>'."\n\n";

			## Counting for required colors
			my $ref_sequence_count = scalar (keys %{$sequences{$reference}});
			if ($orientation eq 'normal'){
				print "\n"."Total # of sequences in reference $reference = $ref_sequence_count"."\n";
			}

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
			## not just the start
			my $increment = scalar(@color_set)/$ref_sequence_count;
			my $rounded_increment = round($increment);
			my $color_start = 0;

			## Making sure that the increment is >= 1
			if ($rounded_increment == 0){
				$rounded_increment = 1;
			}
			
			foreach my $refseq (sort (keys %{$sequences{$reference}})){
				foreach my $queseq (sort (keys %sequences)){
					if (($queseq eq 'concatenated') or ($queseq eq $reference)){
						next;
					}
					else{
						foreach my $queseq_cg (sort (keys %{$sequences{$queseq}})){
							print $cg '<rule>'."\n";
							print $cg "condition  = between($refseq,$queseq_cg)"."\n";
							print $cg "color      = ".$color_set[$color_start]."\n";
							print $cg 'flow       = continue'."\n";
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

			print $cg '</link>'."\n";
			print $cg '</links>'."\n\n";
		}

		## image.conf
		print $cg '<image>'."\n";
		print $cg '<<include etc/image.conf>>'."\n";
		print $cg '</image>'."\n\n";
		print $cg '<<include etc/colors_fonts_patterns.conf>>'."\n";
		print $cg '<<include etc/housekeeping.conf>>'."\n";

		close $cg;
	
	}

}

####################################################################
### Running Circos
####################################################################

## Running circos; sometimes breaks because it doesn't
## find the links file (even if present); runs fine if done manually
## sleep() timer doesn't seem to fix it...
## Not sure what causes Circos to misbehave

if ($circos_plot){

	sleep(10);

	for my $genome ((keys %sequences), 'concatenated'){

		my $subdir = $outdir.'/'.$genome;
		my $subfile = $subdir.'/'.$genome;
		my $config = $subfile.'.conf';
		my $image = $genome.'.png';

		my $pngdir = $subdir.'/images/';
		unless (-d $pngdir){
			make_path($pngdir,{mode => 0755}) or die "Can't create $pngdir: $!\n";
		}

		for my $orientation ('normal', 'inverted'){

			if ($orientation eq 'inverted'){
				$config = $subfile.'.inverted.conf';
				$image = $genome.'.inverted.png';
			}

			print "\n\nRunning Circos for $genome: $orientation\n\n";

			system ("circos \\
			  -conf $config \\
			  -outputfile $image \\
			  -outputdir $pngdir"
			);
		
		}

	}

}

####################################################################
### Subroutine(s)
####################################################################

sub biases {

		my $curseq = $_[0];
		my $pos = $_[1];
		my $divider = $_[2];

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

sub custom_colors {

	my $color_file = $outdir.'/custom_colors.conf';
	open COL, '>', $color_file or die "Can't create $color_file: $!\n";

	print COL '### Custom palette of 20 colors created for the Chloropicon manuscript:'."\n";
	print COL '### https://www.nature.com/articles/s41467-019-12014-x'."\n\n";

	print COL '### Add these colors to Circos etc/colors.conf'."\n";
	print COL 'c01	= 202,75,75'."\n";
	print COL 'c02	= 239,60,104'."\n";
	print COL 'c03	= 241,102,140'."\n";
	print COL 'c04	= 245,152,162'."\n";
	print COL 'c05	= 245,126,47'."\n";
	print COL 'c06	= 250,166,55'."\n";
	print COL 'c07	= 255,197,61'."\n";
	print COL 'c08	= 255,228,67'."\n";
	print COL 'c09	= 210,213,76'."\n";
	print COL 'c10	= 147,195,84'."\n";
	print COL 'c11	= 18,178,89'."\n";
	print COL 'c12	= 0,179,127'."\n";
	print COL 'c13	= 0,180,161'."\n";
	print COL 'c14	= 0,182,204'."\n";
	print COL 'c15	= 0,183,241'."\n";
	print COL 'c16	= 0,157,218'."\n";
	print COL 'c17	= 64,131,196'."\n";
	print COL 'c18	= 94,104,176'."\n";
	print COL 'c19	= 108,82,162'."\n";
	print COL 'c20	= 122,42,144'."\n";

}