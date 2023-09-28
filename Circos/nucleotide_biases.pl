#!/usr/bin/perl
## Pombert Lab, 2022
my $name = 'nucleotide_biases.pl';
my $version = '0.4';
my $updated = '2023-09-28';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;
use File::Path qw(make_path);

my $usage = <<"OPTIONS";
NAME		$name
VERSION		$version
UPDATED		$updated
SYNOPSIS	Generates tab-delimited sliding windows of GC, AT, purine and pyrimidine
		 distributions for easy plotting with MS Excel or other tool. Can also generate
		 coordinate files for Circos.

COMMAND		$name \\
		  -f *.fasta \\
		  -o output_directory \\
		  -w 1000 \\
		  -s 500 \\
		  -c \\
		  -g 0 \\
		  -i \\
		  -r CCMP1205

-f (--fasta)	Fasta file(s) to process
-o (--outdir)	Output directory [Default: ntBiases]
-w (--winsize)	Sliding window size [Default: 10000]
-s (--step)		Sliding window step [Default: 5000]
-c (--circos)	Output files for Circos plotting
-r (--reference)	Genome reference for Circos plotting
-g (--gap)		Default gap links file for Circos plotting [Default: 0]
-i (--image)	Generate Circos images
-t (--tsv)		Output tab-delimited files (e.g. for excel plotting)
OPTIONS
die "\n$usage\n" unless @ARGV;

my @fasta;
my $outdir = 'ntBiases';
my $winsize = 1000;
my $step = 500;
my $circos;
my $reference;
my $gap = 0;
my $image;
my $tsv;
GetOptions(
	'f|fasta=s@{1,}' => \@fasta,
	'o|outdir=s' => \$outdir,
	'w|winsize=i' => \$winsize,
	's|step=i' => \$step,
	'c|circos' => \$circos,
	'r|reference=s' => \$reference,
	'g|gap=i' => \$gap,
	'i|image' => \$image,
	't|tsv' => \$tsv
);

### Check if output directory / subdirs can be created
$outdir =~ s/\/$//;
unless (-d $outdir) {
	make_path($outdir,{mode => 0755})  or die "Can't create $outdir: $!\n";
}
my $catdir = $outdir.'/CONCATENATED';
unless (-d $catdir) {
	make_path($catdir,{mode => 0755}) or die "Can't create $catdir: $!\n";
}

### Hash to store percentage values
my %percent;

### Creating concatenated files
my $cat_kar = $catdir.'/concatenated.genotype';
open CONCAT, ">", $cat_kar or die "Can't create $cat_kar: $!\n";
print CONCAT '#chr - ID LABEL START END COLOR'."\n";

no warnings 'once';
my @cathandles = (*CGC, *CAT, *CGT, *CAC, *CGA, *CCT);
foreach my $cfh (@cathandles){
	my ($lfh) = $cfh =~ /C(\w+)$/;
	my $cat_bias = $catdir.'/concatenated.'.$lfh;
	open $cfh, '>', $cat_bias or die "Can't create $cat_bias: $!\n";
	print $cfh '#chr START END GC_PERCENTAGE'."\n";
}

### Iterating through FASTA file(s)

my %data;
$data{'CONCATENATED'} = 1;

### Creating database of sequences (could be multifasta)
my %sequences;
my $seqname;
my $fileprefix;

while (my $fasta = shift@fasta){

	my ($basename) = fileparse($fasta);
	($fileprefix) = $basename =~ /(\S+)\.\w+$/;

	$data{$fileprefix} = 1;

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

	my $circos_GC;
	my $circos_AT;
	my $circos_GT;
	my $circos_AC;
	my $circos_GA;
	my $circos_CT;
	
	my $circos_kar;
	my $id = 0;

	if ($circos){

		my $circos_dir = $outdir.'/'.$fileprefix;
		unless (-d $circos_dir) {
				make_path($circos_dir,{mode => 0755}) or die "Can't create $circos_dir: $!\n";
		}

		## Creating a "karyotype" file for Circos
		$circos_kar = $circos_dir.'/'.$fileprefix.'.genotype';
		open KAR, ">", $circos_kar or die "Can't create $circos_kar: $!\n";
		print KAR '#chr - ID LABEL START END COLOR'."\n";
		
		## Nucleotide biases
		for my $fh (@filehandles){
			my ($lfh) = $fh =~ /(\w+)$/; ## grabbing test from filehandle
			my $filename = $circos_dir.'/'.$fileprefix.'.'.$lfh;
			open $fh, ">", $filename or die "Can't create $filename: $!\n";
			print $fh '#chr START END GC_PERCENTAGE'."\n";
		}

	}

	### Iterating through each sequence in the FASTA file
	foreach my $sequence (sort (keys %{$sequences{$fileprefix}})){

		my $outfile = $fasta_dir.'/'.$fileprefix.'.'.$sequence.'.tsv';

		if ($tsv){
			open BIAS, ">", $outfile or die "Can't create $outfile: $!\n";
			print BIAS "# Location\t% GC\t% AT\t% AG\t% CT\t% GT\t% AC\n";
		}

		### Sliding windows
		my $seq = $sequences{$fileprefix}{$sequence};
		my $csize = length $seq;

		if ($circos){
			my $terminus = $csize - 1;
			$id++;
			print KAR "chr - $sequence $id 0 $terminus black\n";
			print CONCAT "chr - $sequence $id 0 $terminus black\n";
		}

		my $x;
		for ($x = 0; $x <= ($csize - $winsize); $x += $step){

			my $subseq = substr($seq, $x, $winsize);
			my $end = $x + $winsize - 1;
			
			%percent = ();
			biases($subseq,$x);

			if ($circos){
				foreach my $fh (@filehandles){
					my ($lfh) = $fh =~ /(\w+)$/;
					print $fh "$sequence $x $end $percent{$lfh}\n";
				}
				foreach my $cfh (@cathandles){
					my ($lfh) = $cfh =~ /C(\w+)$/;
					print $cfh "$sequence $x $end $percent{$lfh}\n";
				}
			}

		}

		### Working on leftover string < $winsize
		my $modulo = $csize % $winsize;
		my $subseqleft = substr($seq, -$modulo, $modulo);
		my $leftover_size = length $subseqleft;

		%percent = ();
		biases($subseqleft,$x);

		if ($circos){
			my $end = $csize - 1;
			foreach my $fh (@filehandles){
				my ($lfh) = $fh =~ /(\w+)$/;
				print $fh "$sequence $x $end $percent{$lfh}\n";
			}
			foreach my $cfh (@cathandles){
				my ($lfh) = $cfh =~ /C(\w+)$/;
				print $cfh "$sequence $x $end $percent{$lfh}\n";
			}
		}

		close BIAS;

	}

	close FASTA;

}

####################################################################
### Circos
####################################################################

## Ideogram for Circos
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

############## Ticks
my $ticks = $outdir.'/'.'ticks.conf';
open my $tk, '>', $ticks or die $!;

my $ticks_data = <<'TICKS';
show_ticks          = yes
show_tick_labels    = yes

<ticks>
radius           = 1r
orientation	= out
tick_separation  = 3p
label_separation = 1p
multiplier       = 1e-6
color            = black
thickness        = 3p

<tick>
size             = 20p
spacing        = 10u
show_label     = yes
label_size     = 15p
label_offset   = 5p
format         = %d
suffix         = kb
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
format         = %d
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

my $axes = <<'AXES';
<axes>
	<axis>
	color     = lgrey_a2
	thickness = 1
	spacing   = 0.025r
	</axis>
</axes>
AXES

my @biases = ('GC','AT','GT','AC','GA','CT');
my %colors = (
	'GC' => 'vdred',
	'AT' => 'vdgrey',
	'GT' => 'vdblue',
	'AC' => 'vdgreen',
	'GA' => 'vdpurple',
	'CT' => 'vdyellow'
);

############## conf
for my $genome (keys %data){

	my $subdir = $outdir.'/'.$genome;

	if ($genome eq 'CONCATENATED'){
		$genome = lc($genome);
	}

	my $karyotype = $subdir.'/'.$genome.'.genotype';
	my $config = $subdir.'/'.$genome.'.conf';
	my $subfile = $subdir.'/'.$genome;
	
	my $image = $genome.'.png';

	## Creating a conf file for Circos
	open my $cg, '>', $config or die $!;

	print $cg "<<include $ideogram>>"."\n";
	print $cg "<<include $ticks>>"."\n\n";
	print $cg "karyotype = $karyotype"."\n";
	print $cg 'chromosomes_units=10000'."\n\n";

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

		print $cg $axes;
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
		my $color_counter = sprintf("%02d",1);
		foreach my $refseq (sort (keys %{$sequences{$reference}})){
			foreach my $queseq (sort (keys %sequences)){
				if (($queseq eq 'concatenated') or ($queseq eq $reference)){
					next;
				}
				else{
					foreach my $queseq_cg (sort (keys %{$sequences{$queseq}})){
						print $cg '<rule>'."\n";
						print $cg "condition  = between($refseq,$queseq_cg)"."\n";
						print $cg "color      = c$color_counter"."\n";
						print $cg 'flow       = continue'."\n";
						print $cg '</rule>'."\n\n";
					}
				}
			}
			$color_counter++;
			$color_counter = sprintf("%02d",$color_counter);
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

	## Running circos; currently breaks due to naming
	## Scheme for the links file; need to rethink it

	if ($image){

		my $pngdir = $subdir.'/images/';
		unless (-d $pngdir){
			make_path($pngdir,{mode => 0755}) or die "Can't create $pngdir: $!\n";
		}

		print "\n\nRunning Circos on $pngdir\n\n";
		system "circos \\
		-conf $config \\
		-outputfile $image \\
		-outputdir $pngdir";
	
	}

}

### colors
my $custom_colors =<<'COLORS';
c01	= 202,75,75
c02	= 239,60,104
c03	= 241,102,140
c04	= 245,152,162
c05	= 245,126,47
c06	= 250,166,55
c07	= 255,197,61
c08	= 255,228,67
c09	= 210,213,76
c10	= 147,195,84
c11	= 18,178,89
c12	= 0,179,127
c13	= 0,180,161
c14	= 0,182,204
c15	= 0,183,241
c16	= 0,157,218
c17	= 64,131,196
c18	= 94,104,176
c19	= 108,82,162
c20	= 122,42,144

hm-1 = 245,255,255
hm-2 = 190,233,244
hm-3 = 180,212,233
hm-4 = 170,191,222
hm-5 = 160,171,211
hm-6 = 149,150,200
hm-7 = 138,130,190
hm-8 = 127,111,179
hm-9 = 115,91,168
hm-10 = 103,71,157
hm-11 = 91,52,146
hm-12 = 77,30,136
hm-13 = 63,0,125
COLORS

my $color_file = $outdir.'/custom_colors.conf';
open COL, '>', $color_file or die "Can't create $color_file: $!\n";
print COL $custom_colors;

### Subroutine(s)

sub biases {

		my $curseq = $_[0];
		my $pos = $_[1];

		my $gc = $curseq =~ tr/GgCc//;
		my $at = $curseq =~ tr/AaTt//;
		my $ga = $curseq =~ tr/GgAa//;
		my $ct = $curseq =~ tr/CcTt//;
		my $gt = $curseq =~ tr/GgTt//;
		my $ac = $curseq =~ tr/AaCc//;
		

		$percent{'GC'} = $gc = ($gc/$winsize) * 100;
		$percent{'AT'} = $at = ($at/$winsize) * 100;
		$percent{'GA'} = $ga = ($ga/$winsize) * 100;
		$percent{'CT'} = $ct = ($ct/$winsize) * 100;
		$percent{'GT'} = $gt = ($gt/$winsize) * 100;
		$percent{'AC'} = $ac = ($ac/$winsize) * 100;

		if ($tsv){
			print BIAS "$pos\t$gc\t$at\t$ga\t$ct\t$gt\t$ac\n";
		}

}