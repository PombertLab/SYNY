#!/usr/bin/env perl
# Pombert lab, 2022

my $name = 'run_syny.pl';
my $version = '0.7g';
my $updated = '2024-05-22';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Getopt::ArgvFile home=>1;
use File::Basename;
use Cwd qw(abs_path);
use File::Path qw(make_path);
use threads;
use threads::shared;

###################################################################################################
## Options
###################################################################################################

my $usage = <<"EXIT";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Runs the SYNY pipeline

REQS		Perl / PerlIO::gzip - https://metacpan.org/pod/PerlIO::gzip
		Python3 / matplotlib / seaborn / pandas - https://www.python.org/
		Diamond - https://github.com/bbuchfink/diamond
		Minimap2 - https://github.com/lh3/minimap2
		Circos - http://circos.ca/

USAGE		${name} \\
		  -t 16 \\
		  -a *.gbff \\
		  -o SYNY \\
		  -e 1e-10 \\
		  -g 0 1 5 \\
		  --circos pair

OPTIONS:
-h (--help)             Display all command line options
-t (--threads)          Number of threads to use [Default: 16]
-p (--pthreads)         Number of graphs to plot in parralel; defaults to --threads if unspecified
-a (--annot)            GenBank GBF/GBFF Annotation files (GZIP files are supported)
-o (--outdir)           Output directory [Default = SYNY]
-e (--evalue)           DIAMOND BLASTP evalue cutoff [Default = 1e-10]
-g (--gaps)             Allowable number of gaps between gene pairs [Default = 0]

--minsize               Minimum contig size (in bp) [Default: 1]
--include               Select contigs with names from input text file(s) (one name per line); i.e. exclude everything else
--ranges                Select contigs with subranges from input text file(s): name start end
--exclude               Exclude contigs with names matching the regular expression(s); e.g. --exclude '^AUX'
--aligner               Specify genome alignment tool: minimap or mashmap [Default: minimap]
--asm                   Specify minimap max divergence preset (--asm 5, 10 or 20) [Default: off]
--mpid                  Specify mashmap3 percentage identity [Default: 85]
--resume                Resume minimap/mashmap computations (skip completed alignments)
--no_map                Skip minimap/mashmap pairwise genome alignments
--no_clus               Skip gene cluster reconstructions
--version               Display SYNY version
EXIT

my $plot_options = <<"PLOT_OPTIONS";
### Circos plots
-c (--circos)           Circos plot mode: pair (pairwise), cat (concatenated), all (cat + pair) [Default: pair]
--orientation           Karyotype orientation: normal, inverted or both [Default: normal]
--circos_prefix         Prefix for concatenated plots [Default: circos]
-r (--ref)              Reference to use for concatenated plots; uses first genome (alphabetically) if ommitted
-u (--unit)             Size unit (Kb or Mb) [Default: Mb]
--winsize               Sliding windows size (nucleotide biases) [Default: 10000]
--stepsize              Sliding windows step (nucleotide biases) [Default: 5000]
--labels                Contig label type: mixed (arabic + roman numbers), arabic, roman, or names [Default: mixed]
--label_size            Contig label size [Default: 36]
--label_font            Contig label font [Default: bold]
--custom_file           Load custom colors from file
--list_preset           List available custom color presets
--custom_preset         Use a custom color preset, e.g.: --custom_preset chloropicon
--max_ticks             Set max number of ticks [Default: 5000]
--max_ideograms         Set max number of ideograms [Default: 200]
--max_links             Set max number of links [Default: 75000]
--max_points_per_track  Set max number of points per track [Default: 75000]
--clusters              Color by cluster instead of contig/chromosome [Default: off]
--no_ntbiases           Turn off nucleotide biases subplots
--no_cticks             Turn off ticks in Circos plots
--no_circos             Turn off Circos plots

### Barplots
-bh (--bheight)         Barplot figure height in inches [Default: 10.8]
-bw (--bwidth)          Barplot figure width in inches [Default: 19.2]
--bfsize                Barplot font size [Default: 8]
--palette               Barplot color palette [Default: Spectral]
--monobar               Use a monochrome barplot color instead: e.g. --monobar blue
--no_barplot            Turn off barplots

### Dotplots
-dh (--dheight)         Dotplot figure height in inches [Default: 10.8]
-dw (--dwidth)          Dotplot figure width in inches [Default: 19.2]
--dfsize                Dotplot font size [Default: 8]
-m (--multi)            Axes units multiplier (for dotplots) [Default: 1e5]
--color                 Dotplot color [Default: blue]
--dotpalette            Use a color palette instead: e.g. --dotpalette inferno
--noticks               Turn off ticks on x and y axes
--wdis                  Horizontal distance (width) between subplots [Default: 0.05]
--hdis                  Vertical distance (height) between subplots [Default: 0.1]
--no_dotplot            Turn off dotplots

### Heatmaps
-hh (--hheight)         Heatmap figure height in inches [Default: 10]
-hw (--hwidth)          Heatmap figure width in inches [Default: 10]
--hfsize                Heatmap font size [Default: 8]
--hmpalette             Heatmap color palette [Default: winter_r]
--hmax                  Set maximum color bar value [Default: 100]
--hmin                  Set minimum color bar value [Default: 0]
--hauto                 Set color bar values automatically instead
--no_heatmap            Turn off heatmaps

### Linear maps
-lh (--lheight)         Linear map figure height in inches [Default: 5]
-lw (--lwidth)          Heatmap figure width in inches [Default: 20]
--lm_rpalette           Reference genome color palette [Default: tab20]
--lm_xpalette           Target genome color palette [Default: Blues]
--lmrotation            Contig name rotation [Default: 90]
--lfsize                Font size [Default: 8]
--no_linemap            Turn off linemaps
PLOT_OPTIONS

die ("\n$usage\n") unless (@ARGV);

my @commands = @ARGV;

# Main
my @annot_files;
my $evalue = '1e-10';
my @gaps;
my $outdir = 'SYNY';
my $threads = 16;
my $max_pthreads = $threads;
my $minsize = 1;
my @excluded;
my @included;
my @ranges;
my $aligner = 'minimap';
my $nomap;
my $noclus;
my $resume;
my $asm;
my $mashmap_pid = 85;
my $help;
my $syny_version;

# Circos
my $reference;
my $unit = 'Mb';
my $labels = 'mixed';
my $label_size = 36;
my $label_font = 'bold';
my $circos = 'pair';
my $circos_prefix = 'circos';
my $winsize = 10000;
my $stepsize = 5000;
my $custom_file;
my $custom_colors;
my $list_preset;
my @formats;
my $max_ticks = 5000;
my $max_ideograms = 200;
my $max_links = 75000;
my $max_points_per_track = 75000;
my $clusters;
my $circos_orientation = 'normal';
my $no_ntbiases;
my $no_circos;
my $no_cticks;

# Barplots
my $bheight = 10.8;
my $bwidth = 19.2;
my $bfsize = 8;
my $palette = 'Spectral';
my $monobar;
my $no_barplot;

# Dotplots
my $dheight = 10.8;
my $dwidth = 19.2;
my $dfsize = 8;
my $multiplier = '1e5';
my $color = 'blue';
my $dotpalette;
my $noticks;
my $wdis = 0.05;
my $hdis = 0.1;
my $no_dotplot;

# Heatmaps
my $hheight = 10;
my $hwidth = 10;
my $hfsize = 8;
my $hmpalette = 'winter_r';
my $hmax = 100;
my $hmin = 0;
my $hauto;
my $no_heatmap;

# Linemaps
my $lheight = 5;
my $lwidth = 20;
my $lm_rpalette = 'tab20';
my $lm_xpalette = 'Blues';
my $lmrotation = 90;
my $lfsize = 8;
my $no_linemap;

GetOptions(
	# Main
	't|threads=i' => \$threads,
	'p|pthreads=i' => \$max_pthreads,
	'a|annot=s@{1,}' => \@annot_files,
	'o|outdir=s' => \$outdir,
	'e|evalue=s' => \$evalue,
	'g|gaps=i{0,}' => \@gaps,
	'minsize=i' => \$minsize,
	'excluded=s{0,}' => \@excluded,
	'included=s{0,}' => \@included,
	'ranges=s{0,}' => \@ranges,
	'aligner=s' => \$aligner,
	'no_map' => \$nomap,
	'no_clus' => \$noclus,
	'resume' => \$resume,
	'asm=i' => \$asm,
	'mpid=s' => \$mashmap_pid,
	'h|help' => \$help,
	'version' => \$syny_version,
	# Circos
	'c|circos=s' => \$circos,
	'r|ref|reference=s' => \$reference,
	'u|unit=s' => \$unit,
	'labels=s' => \$labels,
	'label_size=s' => \$label_size,
	'label_font=s' => \$label_font,
	'winsize=i' => \$winsize,
	'stepsize=i' => \$stepsize,
	'circos_prefix=s' => \$circos_prefix,
	'custom_file=s' => \$custom_file,
	'custom_preset=s' => \$custom_colors,
	'list_preset'	=> \$list_preset,
	'f|format=s{1,}' => \@formats,
	'max_ticks=i' => \$max_ticks,
	'max_ideograms=i' => \$max_ideograms,
	'max_links=i' => \$max_links,
	'max_points_per_track=i' => \$max_points_per_track,
	'clusters' => \$clusters,
	'orientation=s' => \$circos_orientation,
	'no_ntbiases' => \$no_ntbiases,
	'no_circos' => \$no_circos,
	'no_cticks' => \$no_cticks,
	# Barplots
	'bh|bheight=s' => \$bheight,
	'bw|bwidth=s' => \$bwidth,
	'bfsize=i' => \$bfsize,
	'palette=s' => \$palette,
	'monobar=s' => \$monobar,
	'no_barplot' => \$no_barplot,
	# Dotplots
	'dh|dheight=s' => \$dheight,
	'dw|dwidth=s' => \$dwidth,
	'dfsize=i' => \$dfsize,
	'm|multiplier=s' => \$multiplier, 
	'color=s' => \$color,
	'dotpalette=s' => \$dotpalette,
	'noticks' => \$noticks,
	'wdis=s' => \$wdis,
	'hdis=s' => \$hdis,
	'no_dotplot' => \$no_dotplot,
	# Heatmaps
	'hh|hheight=s' => \$hheight,
	'hw|hwidth=s' => \$hwidth,
	'hfsize=i' => \$hfsize,
	'hmpalette=s' => \$hmpalette,
	'hmax=i' => \$hmax,
	'hmin=i' => \$hmin,
	'hauto' => \$hauto,
	'no_heatmap' => \$no_heatmap,
	# Linemaps
	'lh|lheight=s' => \$lheight,
	'lw|lwidth=s' => \$lwidth,
	'lm_rpalette=s' => \$lm_rpalette,
	'lm_xpalette=s' => \$lm_xpalette,
	'lmrotation=i' => \$lmrotation,
	'lfsize=i' => \$lfsize,
	'no_linemap' => \$no_linemap
);

# Displaying the full list of options
if ($help){
	print "\n".$usage."\n";
	print $plot_options."\n";
	exit;
}

# Setting default gap value
unless(@gaps){
	@gaps = (0);
}

# Grabbing $path location from script
my ($script,$path) = fileparse($0);
my $align_path = $path.'/Alignments';
my $cluster_path = $path.'/Clusters';
my $plot_path = $path.'/Plots';
my $util_path = $path.'/Utils';

# plot threads
my $pthreads = $threads;
if ($max_pthreads){
	$pthreads = $max_pthreads;
}

###################################################################################################
## Display syny version
###################################################################################################

if ($syny_version){
	my $version_file = $path.'./versions.txt';
	open VER, '<', $version_file or die "Can't read $version_file: $!\n";
	print "\n";
	while (my $line = <VER>){
		print $line;
	}
	print "\n\n";
	exit;
}

###################################################################################################
## Precheck Circos options
###################################################################################################

if ($list_preset){

	system("
		$plot_path/nucleotide_biases.pl \\
		--list_preset
	") == 0 or checksig();
	exit;

}

## Check for plot orientation
my %circos_orientations = (
	'normal' => '',
	'inverted' => '',
	'both' => ''
);

unless ($no_circos){
	$circos_orientation = lc($circos_orientation);
	if (!exists $circos_orientations{$circos_orientation}){
		print "\nCircos orientation $circos_orientation not found. Possible oritentations are:"."\n\n";
		for my $key (sort (keys %circos_orientations)){
			print '  '.$key."\n";
		}
		print "\n";
		exit;
	}

}

## Check for plot fonts
my %circos_fonts = (
	'light' => '',
	'normal' => '',
	'default' => '',
	'semibold' => '',
	'bold' => '',
	'italic' => '',
	'bolditalic' => '',
	'italicbold' => ''
);

if ($label_font){
	my $font = lc($label_font);
	if (!exists $circos_fonts{$font}){
		print "\nCircos font $font not found. Possible fonts are:"."\n\n";
		for my $key (sort (keys %circos_fonts)){
			print '  '.$key."\n";
		}
		print "\n";
		exit;
	}
}

###################################################################################################
## Precheck matplotlib colors
###################################################################################################

my @matp_colors = ($palette, $hmpalette, $lm_rpalette, $lm_xpalette);

if ($dotpalette){
	push (@matp_colors, $dotpalette);
}

my $matp_check = `$util_path/check_mp_colors.py --check @matp_colors`;

if ($matp_check =~ /\[E\]/){
	print $matp_check;
	exit();
}

###################################################################################################
## Checking dependencies
###################################################################################################

# Diamond
my $diamond_check = `echo \$(command -v diamond)`;
chomp $diamond_check;
if ($diamond_check eq ''){
	print STDERR "\n[E]: Cannot find diamond. Please install diamond in your \$PATH. Exiting..\n\n";
	exit;
}

# minimap2
if ($aligner =~ /minimap/i){
	my $minimap2_check = `echo \$(command -v minimap2)`;
	chomp $minimap2_check;
	if ($minimap2_check eq ''){
		print STDERR "\n[E]: Cannot find minimap2. Please install minimap2 in your \$PATH. Exiting..\n\n";
		exit;
	}
}

# mashmap3
if ($aligner =~ /mashmap/i){
	my $mashmap_check = `echo \$(command -v mashmap)`;
	chomp $mashmap_check;
	if ($mashmap_check eq ''){
		print STDERR "\n[E]: Cannot find mashmap. Please install mashmap in your \$PATH. Exiting..\n\n";
		exit;
	}
}

# Circos
unless ($no_circos){
	my $circos_check = `echo \$(command -v circos)`;
	chomp $circos_check;
	if ($diamond_check eq ''){
		print STDERR "\n[E]: Cannot find circos. Please install circos in your \$PATH. Exiting..\n\n";
		exit;
	}
}

###################################################################################################
## Precheck annotations files
###################################################################################################

if (!@annot_files){
	print "\n[E] No annotation file entered with -a/--annot.\n";
	print "[E] Please enter some annotation files to compare. Exiting...\n\n";
	exit;
}
else {

	my $gbff_counter = 0;
	my @files;

	for my $file (@annot_files){
		if ($file =~ /\.gbff(\.gz)?$/){
			$gbff_counter++;
			push (@files, $file);
		}
	}

	if ($gbff_counter == 0){
		print "\n[E] No .gbff/.gbff.gz file detected: --annot @annot_files\n";
		print "[E] Please check the command line. Exiting...\n\n";
		exit;
	}
 	elsif ($gbff_counter == 1){
		print "\n[E] Only one annotation file (.gbff/.gbff.gz) detected: $files[0]\n";
		print "[E] Please enter at least two annotation files to compare. Exiting...\n\n";
		exit;
	}
}

###################################################################################################
## Output directory creation and setup
###################################################################################################

## Creating path before running abs_path(); otherwise returns nothing
unless (-d $outdir){
	make_path($outdir,{mode=>0755}) or die "Can't create $outdir: $!\n";
}
$outdir = abs_path($outdir);

my $mmap_dir = $outdir.'/ALIGNMENTS';
my $list_dir = $outdir.'/LISTS';

## Sequence subdirs
my $seq_dir = $outdir.'/SEQUENCES';
my $prot_dir = $seq_dir.'/PROTEINS';
my $genome_dir = $seq_dir.'/GENOMES';

# Gene cluster subdirs
my $cluster_dir = $outdir.'/CLUSTERS';
my $conserved_dir = $cluster_dir.'/HOMOLOGS';
my $diamond_dir = $cluster_dir.'/DIAMOND';
my $db_dir = $diamond_dir.'/DB';
my $cluster_synteny = $cluster_dir.'/SYNTENY';

# Plots subdirs
my $plots_dir = $outdir.'/PLOTS';
my $paf_hm_dir = $plots_dir.'/HEATMAPS';
my $barplot_dir = $plots_dir.'/BARPLOTS';
my $dotplot_dir = $plots_dir.'/DOTPLOTS';
my $linemap_dir = $plots_dir.'/LINEMAPS';
my $circos_plot_dir = $plots_dir.'/CIRCOS';
my $circos_data_dir = $plots_dir.'/CIRCOS_DATA';
my $circos_cat_dir = $circos_data_dir.'/concatenated';
my $circos_pair_dir = $circos_data_dir.'/pairwise';

my @outdirs = (
	$list_dir,
	$seq_dir,$prot_dir,
	$genome_dir,
	$plots_dir,
	$circos_data_dir,
	$circos_cat_dir,
	$circos_pair_dir,
	$diamond_dir,
	$db_dir
);

unless ($noclus){
	push (@outdirs, $cluster_dir);
	push (@outdirs, $cluster_synteny);
}

foreach my $dir (@outdirs){
	unless (-d $dir){
		make_path($dir,{mode=>0755}) or die "Can't create $dir: $!\n";
	}
}

## Logs
my $log_file = $outdir.'/'.'syny.log';
my $log_err = $outdir.'/'.'error.log';

open LOG, '>', $log_file or die "Can't create $log_err: $!\n";
open ERROR, '>', $log_err or die "Can't create $log_err: $!\n";

my $time = localtime();
my $start_time = time();
my $tstart = time();

print LOG "SYNY started on: ".$time."\n\n";
logs(\*LOG, 'Command line', 'header');
print LOG "$0 @commands\n\n";

###################################################################################################
## Options flags
###################################################################################################

my $excluded_flag = '';
if (@excluded){
	$excluded_flag = "--exclude @excluded";
}

my $included_flag = '';
if (@included){
	$included_flag = "--include @included";
}

my $ranges_flag = '';
if (@ranges){
	$ranges_flag = "--ranges @ranges";
}


my $cluster_flag = '';
if ($clusters){
	$cluster_flag = '--clusters';
}

my $custom_cc_file = '';
if ($custom_file){
	$custom_cc_file = "--custom_file $custom_file";
}

my $custom_cc = '';
if ($custom_colors){
	$custom_cc = "--custom_preset $custom_colors";
}

my $tick_flag = '';
if ($noticks){
	$tick_flag = '--noticks'
}

my $monobar_flag = '';
if ($monobar){
	$monobar_flag = "--mono $monobar";
}

my $dotpal_flag = '';
if ($dotpalette){
	$dotpal_flag = "--palette $dotpalette";
}

my $asm_flag = '';
if ($asm){
	$asm_flag = "--asm $asm";
}

my $resume_flag = '';
if ($resume){
	$resume_flag = '--resume';
}

my $hm_vauto_flag = '';
if ($hauto){
	$hm_vauto_flag = '--vauto';
}

###################################################################################################
## Run list_maker.pl
###################################################################################################

my @threads = initThreads($threads);

print ERROR "\n### list_maker.pl ###\n";
print "\n##### Extracting data from GenBank files\n";

logs(\*LOG, 'Parsing data', 'header');
logs(\*LOG, 'list_maker.pl', 'start');

my @annotations :shared = @annot_files;
my $annot_num = scalar(@annot_files);

for my $thread (@threads){
	$thread = threads->create(\&list_maker);
}
for my $thread (@threads){
	$thread ->join();
}

logs(\*LOG, 'list_maker.pl', 'end');

###################################################################################################
## Get PAF files with minimap2
###################################################################################################

# Skip minimap if requested 
if ($nomap){
	goto HOMOLOGY;
}

## Running minimap2/mashmap3 with get_paf.pl
$tstart = time();
print "\n##### Infering colinearity from pairwise genome alignments\n";
print ERROR "\n### get_paf.pl ###\n";

logs(\*LOG, 'Genome alignments', 'header');
logs(\*LOG, 'get_paf.pl', 'start');

my $paf_dir = "$mmap_dir/PAF";

system("
	$align_path/get_paf.pl \\
	  --fasta $genome_dir/*.fasta \\
	  --outdir $mmap_dir \\
	  --aligner $aligner \\
	  --threads $threads \\
	  --percent $mashmap_pid \\
	  $resume_flag \\
	  $asm_flag
") == 0 or checksig();

opendir (PAFDIR, $paf_dir) or die "\n\n[ERROR]\tCan't open $paf_dir: $!\n\n";
while (my $file = readdir(PAFDIR)){
	if ($file =~ /\.paf$/){
		my $paf_file = "$paf_dir/$file";
		if (-z $paf_file){
			print LOG "  - $paf_file is empty. No collinearity or process terminated by OOM killer?\n";
		}
	}
}
closedir PAFDIR;

logs(\*LOG, 'get_paf.pl', 'end');

## Creating Circos links file
$tstart = time();
print ERROR "\n### paf2links.pl ###\n";

my $paf_links = "$circos_cat_dir/concatenated.mmap.links";

logs(\*LOG, 'paf2links.pl', 'start');

system("
	$align_path/paf2links.pl \\
	  --paf $paf_dir \\
	  --links $paf_links \\
	  $cluster_flag \\
	  $custom_cc \\
	  $custom_cc_file \\
	  2>> $log_err
") == 0 or checksig();

logs(\*LOG, 'paf2links.pl', 'end');

#### Calculate/plot PAF metrics
if ($aligner =~ /minimap/){

	$tstart = time();
	print ERROR "\n### paf_metrics.py ###\n";
	print "\n# Calculating metrics (genome alignments):\n";

	my $aln_length_dir = $mmap_dir.'/METRICS';
	my @paf_files;

	opendir (PAFDIR, $paf_dir) or die "\n\n[ERROR]\tCan't open $paf_dir: $!\n\n";
	while (my $file = readdir(PAFDIR)){
		if ($file =~ /\.paf$/){
			push (@paf_files, "$paf_dir/$file");
		}
	}
	closedir PAFDIR;

	logs(\*LOG, 'paf_metrics.py', 'start');

	system ("
		$align_path/paf_metrics.py \\
		--paf @paf_files \\
		--outdir $aln_length_dir \\
		--threads $threads \\
		--height 10.8 \\
		--width 19.2 \\
		--color steelblue \\
		2>> $log_err
	") == 0 or checksig();

	logs(\*LOG, 'paf_metrics.py', 'end');

}

###################################################################################################
## Creating barplots / dotplots / linemaps / heatmaps from minimap2 PAF files
###################################################################################################

# Barplots
unless ($no_barplot){

	$tstart = time();
	print ERROR "\n### paf_to_barplot.py (genome alignments) ###\n";
	print "\n# Barplots (genome alignments):\n";

	logs(\*LOG, 'paf_to_barplot.py', 'start');

	system("
		$plot_path/paf_to_barplot.py \\
		--paf $paf_dir/*.paf \\
		--fasta $genome_dir/*.fasta \\
		--outdir $barplot_dir \\
		--threads $pthreads \\
		--height $bheight \\
		--width $bwidth \\
		--palette $palette \\
		--fontsize $bfsize \\
		$tick_flag \\
		$monobar_flag \\
		2>> $log_err
	") == 0 or checksig();

	logs(\*LOG, 'paf_to_barplot.py', 'end');

}

# Doplots
unless ($no_dotplot){

	$tstart = time();
	print ERROR "\n### paf_to_dotplot.py (genome alignments) ###\n";
	print "\n# Dotplots (genome alignments):\n";

	logs(\*LOG, 'paf_to_dotplot.py', 'start');

	system("
		$plot_path/paf_to_dotplot.py \\
		--paf $paf_dir/*.paf \\
		--fasta $genome_dir/*.fasta \\
		--threads $pthreads \\
		--outdir $dotplot_dir \\
		--unit $multiplier \\
		--height $dheight \\
		--width $dwidth \\
		--color $color \\
		--fontsize $dfsize \\
		$tick_flag \\
		$dotpal_flag \\
		--wdis $wdis \\
		--hdis $hdis \\
		2>> $log_err
	") == 0 or checksig();

	logs(\*LOG, 'paf_to_dotplot.py', 'end');

}

# Linemaps
unless ($no_linemap){

	$tstart = time();
	print ERROR "\n### linear_maps.py (genome alignments) ###\n";
	print "\n# Linemaps (genome alignments):\n";

	logs(\*LOG, 'linear_maps.py', 'start');

	system("
		$plot_path/linear_maps.py \\
		--paf $paf_dir/*.paf \\
		--fasta $genome_dir/*.fasta \\
		--outdir $linemap_dir \\
		--threads $pthreads \\
		--height $lheight \\
		--width $lwidth \\
		--rpalette $lm_rpalette \\
		--xpalette $lm_xpalette \\
		--fontsize $lfsize \\
		--rotation $lmrotation \\
		2>> $log_err
	") == 0 or checksig();

	logs(\*LOG, 'linear_maps.py', 'end');

}

# Heatmaps
unless ($no_heatmap){

	$tstart = time();
	print ERROR "\n### paf_to_hm.py (genome alignments) ###\n";
	print "\n# Heatmaps (genome alignments):\n";

	logs(\*LOG, 'paf_to_hm.py', 'start');

	system("
		$plot_path/paf_to_heatmap.py \\
		--paf $paf_dir/*.paf \\
		--fasta $genome_dir/*.fasta \\
		--outdir $paf_hm_dir \\
		--height $hheight \\
		--width $hwidth \\
		--palette $hmpalette \\
		--matrix $mmap_dir/paf_matrix.tsv \\
		--fontsize $hfsize \\
		--vmax $hmax \\
		--vmin $hmin \\
		$hm_vauto_flag \\
		2>> $log_err
	") == 0 or checksig();

	logs(\*LOG, 'paf_to_hm.py', 'end');

}

###################################################################################################
## Run get_homology.pl
###################################################################################################

# Skip gene cluster if requested
if ($noclus){
	goto CIRCOS;
}

# Else, search for gene clusters
HOMOLOGY:

$tstart = time();
print "\n##### Infering colinearity from protein-coding gene clusters\n\n";
print ERROR "\n### get_homology.pl ###\n";

my @prot_files = ();

opendir (FAA, $prot_dir) or die "\n\n[ERROR]\tCan't open $prot_dir: $!\n\n";
while (my $file = readdir(FAA)){

	if ($file =~ /\.faa$/){

		my $prot_file = $prot_dir.'/'.$file;

		## Skipping blank files
		if (-z $prot_file){
			print "Protein file $prot_file is blank\n";
			next;
		}
		else{
			push (@prot_files, $prot_file);
		}
	}

}
closedir FAA;

## Checking if all protein files are blank, if so no point in searching for shared clusters
if ((scalar @prot_files) == 0){
	if ($noclus){
		print "\n[W] All protein files are blank. No annotation detected.\n\n";
	}
	else{
		print "\n[E] All protein files are blank. Skipping cluster detection...\n\n";
	}
	$noclus = 1;
	goto CIRCOS;
}

## If only one is not blank, skip as well (no point in searching for shared clusters)
if ((scalar @prot_files) == 1){
	print "\n[E] Only one protein file contain data. Skipping cluster detection...\n\n";
	$noclus = 1;
	goto CIRCOS;
}

## Otherwize, proceed
logs(\*LOG, 'Gene clusters', 'header');
logs(\*LOG, 'get_homology.pl', 'start');

system("
	$cluster_path/get_homology.pl \\
	--input @prot_files \\
	--evalue $evalue \\
	--threads $threads \\
	--outdir $diamond_dir \\
	--shared $conserved_dir \\
	--list $list_dir \\
	2>> $log_err
") == 0 or checksig();

logs(\*LOG, 'get_homology.pl', 'end');

###################################################################################################
## Run get_synteny.pl
###################################################################################################

$tstart = time();
print ERROR "\n### get_synteny.pl ###\n";

logs(\*LOG, 'get_synteny.pl', 'start');

my $clu_sum_file = $cluster_synteny.'/'.'clusters_summary.tsv';
if (-f $clu_sum_file){
	system("rm $clu_sum_file") == 0 or checksig();
}

foreach my $prot_file_1 (sort(@prot_files)){

	my ($file_name_1, $dir, $ext) = fileparse($prot_file_1, '\..*');
	print "\nIdentifying synteny between $file_name_1\n"; 

	foreach my $prot_file_2 (sort(@prot_files)){

		if($prot_file_1 ne $prot_file_2){

			my ($file_name_2, $dir, $ext) = fileparse($prot_file_2, '\..*');
			print "\t$file_name_2\n";

			foreach my $gap (@gaps){
				system("
					$cluster_path/get_synteny.pl \\
					--query_list $list_dir/$file_name_1.list \\
					--query_blast $diamond_dir/${file_name_1}_vs_${file_name_2}.diamond.6 \\
					--subject_list $list_dir/$file_name_2.list \\
					--subject_blast $diamond_dir/${file_name_2}_vs_${file_name_1}.diamond.6 \\
					--gap $gap \\
					--outdir $cluster_synteny/gap_$gap \\
					--sumdir $cluster_synteny \\
					2>> $log_err
				") == 0 or checksig();
			}
		}
		else{
			print "\t. . .\n";
		}
	}

}

logs(\*LOG, 'get_synteny.pl', 'end');

###################################################################################################
## Grabbing total number of entries in .list files
###################################################################################################

my %gene_count;
opendir (LISTS, $list_dir) or die "\n\n[ERROR]\tCan't open $list_dir: $!\n\n";

while (my $file = readdir(LISTS)){

	if ($file =~ /\.list$/){

		my $list_file = $list_dir.'/'.$file;
		open TMP, '<', $list_file or die "$!\n";
		my ($basename,$path) = fileparse($file);
		$basename =~ s/\.list$//;

		my $counter = 0;
		while (my $line = <TMP>){
			chomp $line;
			if ($line ne ''){
				$counter++;
			}
		}
		$gene_count{$basename} = $counter;

	}
}

closedir LISTS;

###################################################################################################
## Creating barplots / dotplots / linemaps from SYNY PAF files
###################################################################################################

### Create pseudo PAF files from clusters found
$tstart = time();
print ERROR "\n### clusters_to_paf.pl ###\n";

logs(\*LOG, 'clusters_to_paf.pl', 'start');

foreach my $gap (@gaps){

	my $clusdir = $cluster_synteny.'/gap_'.$gap.'/CLUSTERS';
	my $ppafdir = $cluster_synteny.'/gap_'.$gap.'/PAF';

	system ("
		$cluster_path/clusters_to_paf.pl \\
		--fasta $genome_dir/*.fasta \\
		--lists $list_dir/*.list \\
		--clusters $clusdir/*.clusters \\
		--outdir $ppafdir \\
		2>> $log_err
	") == 0 or checksig();
}

logs(\*LOG, 'clusters_to_paf.pl', 'end');

### Create barplots from PAF files
unless ($no_barplot){

	$tstart = time();
	print ERROR "\n### paf_to_barplot.py (gene clusters) ###\n";
	print "\n# Barplots (gene clusters):\n";

	my @barplot_files = getpaf_files();

	logs(\*LOG, 'paf_to_barplot.py', 'start');

	system("
		$plot_path/paf_to_barplot.py \\
		--paf @barplot_files \\
		--fasta $genome_dir/*.fasta \\
		--threads $pthreads \\
		--outdir $barplot_dir \\
		--height $bheight \\
		--width $bwidth \\
		--palette $palette \\
		--fontsize $bfsize \\
		$tick_flag \\
		$monobar_flag \\
		2>> $log_err
	") == 0 or checksig();

	logs(\*LOG, 'paf_to_barplot.py', 'end');

}

### Create dotplots from PAF files
unless ($no_dotplot){

	$tstart = time();
	print ERROR "\n### paf_to_dotplot.py (gene clusters) ###\n";
	print "\n# Dotplots (gene clusters):\n";

	my @dotplot_files = getpaf_files();

	logs(\*LOG, 'paf_to_dotplot.py', 'start');

	system("
		$plot_path/paf_to_dotplot.py \\
		--paf @dotplot_files \\
		--fasta $genome_dir/*.fasta \\
		--threads $pthreads \\
		--outdir $dotplot_dir \\
		--unit $multiplier \\
		--height $dheight \\
		--width $dwidth \\
		--color $color \\
		--fontsize $dfsize \\
		$tick_flag \\
		$dotpal_flag \\
		--wdis $wdis \\
		--hdis $hdis \\
		2>> $log_err
	") == 0 or checksig();

	logs(\*LOG, 'paf_to_dotplot.py', 'end');

}

### Create linemaps from PAF files
unless ($no_linemap){

	$tstart = time();
	print ERROR "\n### linear_maps.py (gene clusters) ###\n";
	print "\n# Linemaps (gene clusters):\n";

	my @linemap_files = getpaf_files();

	logs(\*LOG, 'linear_maps.py', 'start');

	system("
		$plot_path/linear_maps.py \\
		--paf @linemap_files \\
		--fasta $genome_dir/*.fasta \\
		--outdir $linemap_dir \\
		--threads $pthreads \\
		--height $lheight \\
		--width $lwidth \\
		--rpalette $lm_rpalette \\
		--xpalette $lm_xpalette \\
		--fontsize $lfsize \\
		--rotation $lmrotation \\
		2>> $log_err
	") == 0 or checksig();

	logs(\*LOG, 'linear_maps.py', 'end');

}

### Create a cluster summary table

print ERROR "\n### Cluster summary table ###\n";

my $clu_sum_table = $cluster_synteny.'/'.'clusters_summary_table.tsv';
open CLU, '<', $clu_sum_file or die "Can't read $clu_sum_file: $!\n";
open TSV, '>', $clu_sum_table or die "Can't create $clu_sum_table: $!\n";

print TSV '### Query'."\t";
print TSV 'Total # proteins'."\t";
print TSV 'Allowed Gaps'."\t";
print TSV 'Total # proteins in clusters'."\t";
print TSV '% of proteins in clusters'."\t";
print TSV '# of clusters'."\t";
print TSV 'Longest'."\t".'Shortest'."\t";
print TSV 'Average'."\t".'Median'."\t";
print TSV 'N50'."\t".'N75'."\t".'N90'."\n";

my %cluster_metrics;
my $clu_query;
my $qgap;

my %queries;
my %matrices;

while (my $line = <CLU>){
	chomp $line;
	if ($line =~ /##### (\w+); Gap = (\d+) #####/){
		$clu_query = $1;
		$qgap = $2;
		my ($que) = $clu_query =~ /^(\S+)_vs_(\S+)$/;
		$queries{$que} = '';
	}
	elsif ($line =~ /Total.*\t(\S+)$/){
		$cluster_metrics{$clu_query}{$qgap}{'total'} = $1;
	}
	elsif ($line =~ /# of clusters.*\t(\S+)$/){
		$cluster_metrics{$clu_query}{$qgap}{'total_clu'} = $1;
	}
	elsif ($line =~ /Longest.*\t(\S+)$/){
		$cluster_metrics{$clu_query}{$qgap}{'longest'} = $1;
	}
	elsif ($line =~ /Shortest.*\t(\S+)$/){
		$cluster_metrics{$clu_query}{$qgap}{'shortest'} = $1;
	}
	elsif ($line =~ /Average.*\t(\S+)$/){
		$cluster_metrics{$clu_query}{$qgap}{'average'} = $1;
	}
	elsif ($line =~ /Median.*\t(\S+)$/){
		$cluster_metrics{$clu_query}{$qgap}{'median'} = $1;
	}
	elsif ($line =~ /N50.*\t(\S+)$/){
		$cluster_metrics{$clu_query}{$qgap}{'n50'} = $1;
	}
	elsif ($line =~ /N75.*\t(\S+)$/){
		$cluster_metrics{$clu_query}{$qgap}{'n75'} = $1;
	}
	elsif ($line =~ /N90.*\t(\S+)$/){
		$cluster_metrics{$clu_query}{$qgap}{'n90'} = $1;
	}
}

foreach my $query (sort (keys %cluster_metrics)){
	foreach my $gap (sort {$a <=> $b}(keys %{$cluster_metrics{$query}})){

		print TSV $query."\t";

		my ($ref) = $query =~ /^(\w+)\_vs/;
		print TSV $gene_count{$ref}."\t";

		print TSV $gap."\t";

		print TSV $cluster_metrics{$query}{$gap}{'total'}."\t";
		my $percent = ($cluster_metrics{$query}{$gap}{'total'}/$gene_count{$ref})*100;
		$percent = sprintf("%.2f", $percent);
		print TSV $percent."\t";

		my ($que, $sub) = $query =~ /^(\S+)_vs_(\S+)$/;
		$matrices{$gap}{$que}{$sub} = $percent;

		print TSV $cluster_metrics{$query}{$gap}{'total_clu'}."\t";
		print TSV $cluster_metrics{$query}{$gap}{'longest'}."\t";
		print TSV $cluster_metrics{$query}{$gap}{'shortest'}."\t";
		print TSV $cluster_metrics{$query}{$gap}{'average'}."\t";
		print TSV $cluster_metrics{$query}{$gap}{'median'}."\t";
		print TSV $cluster_metrics{$query}{$gap}{'n50'}."\t";
		print TSV $cluster_metrics{$query}{$gap}{'n75'}."\t";
		print TSV $cluster_metrics{$query}{$gap}{'n90'}."\n";

	}
}

close CLU;
close TSV;

### Create data matrices

print ERROR "\n### Cluster data matrices ###\n";

foreach my $gap (keys %matrices){

	my $matrix_file = $cluster_synteny."/gap_$gap"."/matrix_gap_$gap.tsv";
	open MM, '>', $matrix_file or die $!;

	# header
	foreach my $query (sort (keys %queries)){
		print MM "\t".$query;
	}
	print MM "\n";

	foreach my $query (sort (keys %queries)){
		print MM "$query";
		foreach my $subject (sort (keys %queries)){
			if ($query eq $subject){
				print MM "\t".'100';
			}
			elsif (exists $matrices{$gap}{$query}{$subject}){
				print MM "\t".$matrices{$gap}{$query}{$subject};
			}
			else {
				print MM "\t".'0';
			}
		}
		print MM "\n";
	}

	close MM;

}

### Create cluster summary table as heatmap with matplotlib

unless ($no_heatmap){

	$tstart = time();
	print ERROR "\n### protein_cluster_hm.py (gene clusters) ###\n";
	print "\n# Heatmaps (gene clusters):\n";

	my @hm_files;
	foreach my $gap (@gaps){
		my $matrix_file = $cluster_synteny."/gap_$gap"."/matrix_gap_$gap.tsv";
		push (@hm_files, $matrix_file);
	}

	logs(\*LOG, 'protein_cluster_hm.py', 'start');

	system("
		$plot_path/protein_cluster_hm.py \\
		--tsv @hm_files \\
		--outdir $paf_hm_dir \\
		--threads $pthreads \\
		--height $hheight \\
		--width $hwidth \\
		--palette $hmpalette \\
		--fontsize $hfsize \\
		--vmax $hmax \\
		--vmin $hmin \\
		$hm_vauto_flag \\
		2>> $log_err
	") == 0 or checksig();

	logs(\*LOG, 'protein_cluster_hm.py', 'end');

}

###################################################################################################
## Run id_conserved_regions.pl
###################################################################################################

$tstart = time();
print ERROR "\n### id_conserved_regions.pl ###\n";

system("
	$cluster_path/id_conserved_regions.pl \\
	--lists $list_dir \\
	--blasts $diamond_dir \\
	--outdir $conserved_dir \\
	2>> $log_err
") == 0 or checksig();

logs(\*LOG, 'id_conserved_regions.pl', 'end');

###################################################################################################
## Create links files, karyotypes and nucleotide biases for Circos
###################################################################################################

## Circos links files
$tstart = time();
print ERROR "\n### clusters2links.pl ###\n";

logs(\*LOG, 'clusters2links.pl', 'start');

my @cluster_dirs;
opendir (CDIR, $cluster_synteny) or die "Can't open $cluster_synteny: $!\n";
while (my $dname = readdir(CDIR)) {
	if (-d "$cluster_synteny/$dname"){
		unless (($dname eq '.') or ($dname eq '..')){
			push (@cluster_dirs, "$cluster_synteny/$dname");
		}
	}
}
closedir CDIR;

foreach my $cluster_subdir (@cluster_dirs){
	my ($basedir) = fileparse($cluster_subdir);
	system("
		$cluster_path/clusters2links.pl \\
		--cluster $cluster_subdir/CLUSTERS/*.clusters \\
		--list $list_dir/*.list \\
		--outdir $circos_data_dir \\
		$cluster_flag \\
		$custom_cc_file \\
		$custom_cc \\
		2>> $log_err
	") == 0 or checksig();
}

logs(\*LOG, 'clusters2links.pl', 'end');

###### Circos section
CIRCOS:

## Karyotypes and nucleotide biases
$tstart = time();
print ERROR "\n### nucleotide_biases.pl ###\n";

# Option flags
my $ref = '';
if ($reference){
	$ref = "--reference $reference";
}

my $unit_size = '';
if ($unit){
	$unit_size = "--unit $unit";
}

my $circos_ticks_flag = '';
if ($no_cticks){
	$circos_ticks_flag = '--noticks';
}

my $circos_ntbiases_flag = '';
if ($no_ntbiases){
	$circos_ntbiases_flag = '--no_biases';
}

my $gap = $gaps[0];

## Running nucleotide_biases.pl
system("
	$plot_path/nucleotide_biases.pl \\
	--outdir $circos_data_dir \\
	--fasta $genome_dir/*.fasta \\
	--minsize $minsize \\
	--winsize $winsize \\
	--step $stepsize \\
	--gap $gap \\
	--labels $labels \\
	--label_size $label_size \\
	$ref \\
	$unit_size \\
	$custom_cc_file \\
	$custom_cc \\
	--max_ticks $max_ticks \\
	--max_ideograms $max_ideograms \\
	--max_links $max_links \\
	--max_points_per_track $max_points_per_track \\
	$cluster_flag \\
	$circos_ticks_flag \\
	$circos_ntbiases_flag \\
	2>> $log_err
") == 0 or checksig();

## Create circos configuration files per links file
unless ($noclus){
	foreach my $num (@gaps){
		my $gap = 'gap_'.$num;
		circos_conf($gap);
		pairwise_conf($gap)
	}
}

unless ($nomap){
	circos_conf('mmap');
	pairwise_conf('mmap')
}

##### Running Circos
## Populating list of plots to generate
my %circos_todo_list;

unless ($no_circos){

	logs(\*LOG, 'Circos', 'header');
	logs(\*LOG, 'Plotting', 'start');

	print "\n# Circos plots:\n";

	unless (-d $circos_plot_dir){
			mkdir ($circos_plot_dir, 0755) or die "Can't create $circos_plot_dir: \n";
	}

	## Synteny inferred from protein clusters
	unless ($noclus){
		foreach my $num (@gaps){
			my $affix = 'gap_'.$num;
			if (($circos eq 'cat') or ($circos eq 'all')){
				circos_plot($affix, 'cat');
			}
			if (($circos eq 'pair') or ($circos eq 'all')){
				circos_plot($affix, 'pair');
			}
		}
	}

	## Synteny inferred from minimap2 PAF files
	unless ($nomap){
		my $affix = 'mmap';
		if (($circos eq 'cat') or ($circos eq 'all')){
			circos_plot($affix, 'cat');
		}
		if (($circos eq 'pair') or ($circos eq 'all')){
			circos_plot($affix, 'pair');
		}
	}

}

## Running one circos instance per thread
my @circos_files :shared = sort(keys %circos_todo_list);
my $circos_num = scalar(@circos_files);

unless ($no_circos){

	my @cthreads = initThreads($pthreads);

	for my $thread (@cthreads){
		$thread = threads->create(\&run_circos);
	}
	for my $thread (@cthreads){
		$thread ->join();
	}

	# Moving to PNG/SVG subdirs
	if (($circos eq 'cat') or ($circos eq 'all')){
		my $tmpdir = $circos_plot_dir.'/Concatenated';
		system ("mv $tmpdir/*.png $tmpdir/PNG/");
		system ("mv $tmpdir/*.svg $tmpdir/SVG/");
	}

	if (($circos eq 'pair') or ($circos eq 'all')){
		my $tmpdir = $circos_plot_dir.'/Pairwise';
		system ("mv $tmpdir/*.png $tmpdir/PNG/");
		system ("mv $tmpdir/*.svg $tmpdir/SVG/");
	}

	logs(\*LOG, 'Plotting', 'end');

}


###################################################################################################
## Completion
###################################################################################################

my $end_time = localtime();
my $run_time = time - $start_time;

logs(\*LOG, 'Total runtime', 'header');
print LOG "SYNY completed on: ".$end_time."\n";
print LOG "Total runtime: ".$run_time." seconds\n";
close LOG;

###################################################################################################
## Subroutines
###################################################################################################

sub initThreads {

	my $thread_count = $_[0];
	my @initThreads;

	for (my $i = 1; $i <= $thread_count; $i++){
		push(@initThreads,$i);
	}

	return @initThreads;

}

sub checksig {

	my $exit_code = $?;
	my $modulo = $exit_code % 255;

	if ($modulo == 2) {
		print "\nSIGINT detected: Ctrl+C => exiting...\n\n";
		exit(2);
	}
	elsif ($modulo == 131) {
		print "\nSIGTERM detected: Ctrl+\\ => exiting...\n\n";
		exit(131);
	}

}

sub logs {

	my $fh = $_[0];
	my $analysis = $_[1];
	my $type = $_[2];
	my $len = length($analysis);
	my $pad = 60 - $len;
	my $spacer = ' ' x $pad;


	my $mend = localtime();
	if ($type eq 'header'){
		my $barsize = '#' x 96;
		print $fh $barsize."\n";
		print LOG "##### $analysis\n";
		print LOG $barsize."\n\n";
	}
	elsif ($type eq 'start'){
		print $fh "$analysis:".$spacer."started on $mend\n";
	}
	elsif ($type eq 'end'){

		my $run_time = time - $tstart;
		my $tlen = length($run_time);
		my $tpad = $pad + $len - $tlen - 18;
		my $tspacer = ' ' x $tpad;

		print $fh "Runtime: $run_time seconds".$tspacer."completed on $mend\n\n";

	}

}

sub list_maker {

	while (my $annotation = shift@annotations){

		my $x = $annot_num - scalar(@annotations);
		print "$x / $annot_num - Extracting data from $annotation\n";

		system("
			$util_path/list_maker.pl \\
			--input $annotation \\
			--outdir $outdir \\
			--minsize $minsize \\
			$included_flag \\
			$ranges_flag \\
			$excluded_flag \\
			2>> $log_err
		") == 0 or checksig();

		## Check if annotation file is zero
		my $fname = basename($annotation);
		my @file_data = split('\.',$fname);
		my ($fprefix) = $file_data[0] =~ /^(\w+)/;
		my $prot_dir = $outdir.'/SEQUENCES/PROTEINS';
		my $fasta_proteins = "$prot_dir/$fprefix.faa";
		if (-z $fasta_proteins){
			print LOG ("  - $fname does not contain protein annotations.\n")
		}

	}

	threads->exit();

}

sub getpaf_files {

	my @paf_files;

	foreach my $gap (@gaps){

		my $ppafdir = $cluster_synteny.'/gap_'.$gap.'/PAF';
		opendir (PAFDIR, $ppafdir) or die "\n\n[ERROR]\tCan't open $ppafdir: $!\n\n";

		while (my $file = readdir(PAFDIR)){
			if ($file =~ /\.paf$/){
				push (@paf_files, "$ppafdir/$file");
			}
		}

		closedir PAFDIR;

	}

	return @paf_files;

}

## Circos subs
sub circos_conf {

	my @affix = @_;

	foreach my $affix (shift @affix){

		## Creating Circos configuration file(s) with desired link file
		my $paf_links  = $circos_cat_dir.'/concatenated.'.$affix.'.links';

		for my $orientation ('normal', 'inverted'){

			my $circos_conf_f = $circos_cat_dir.'/concatenated.'.$orientation.'.conf';
			my $new_conf = $circos_cat_dir.'/concatenated.'.$affix.'.'.$orientation.'.conf';

				open CONF, '<', $circos_conf_f or die "Can't read $circos_conf_f: $!\n";
				open NEWCONF, '>', $new_conf or die "Can't create $new_conf: $!\n";

				while (my $line = <CONF>){
					chomp $line;
					if ($line =~ /^file.*.links$/){
						print NEWCONF "file          = $paf_links"."\n";
					}
					else {
						print NEWCONF $line."\n";
					}
				}

				close CONF;
				close NEWCONF;

		}

	}

}

sub pairwise_conf {

	my @affix = @_;

	foreach my $affix (shift @affix){

		## Creating Circos configuration file(s) with desired link file
		my $paf_links = $circos_cat_dir.'/concatenated.'.$affix.'.links';

		opendir(PAIRS, $circos_pair_dir) or die "Can't open $circos_pair_dir: $!\n";
		my @pairs;
		while (my $entry = readdir(PAIRS)){
			if ($entry =~ /_vs_/){
				push(@pairs, $entry);
			}
		}
		closedir PAIRS;

		for my $pair (@pairs){

			for my $orientation ('normal', 'inverted'){

				my $subloc = $circos_pair_dir.'/'.$pair.'/'.$pair;
				my $circos_conf_f = $subloc.'.'.$orientation.'.conf';
				my $pair_genotype = $subloc.'.'.$orientation.'.karyotype';
				my $new_conf = $subloc.'.'.$affix.'.'.$orientation.'.conf';

					open CONF, '<', $circos_conf_f or die "Can't read $circos_conf_f: $!\n";
					open NEWCONF, '>', $new_conf or die "Can't create $new_conf: $!\n";

					while (my $line = <CONF>){
						chomp $line;
						if ($line =~ /^karyotype/){
							print NEWCONF "karyotype = $pair_genotype"."\n";
						}
						elsif ($line =~ /^file.*.links$/){
							print NEWCONF "file          = $paf_links"."\n";
						}
						else {
							print NEWCONF $line."\n";
						}
					}

					close CONF;
					close NEWCONF;

			}

		}

	}

}

sub circos_plot {

	my $affix = $_[0];
	my $mode = $_[1];

	if ($mode eq 'cat'){

		my $plot_dir = $circos_plot_dir.'/Concatenated';
		my $plot_dir_png = $plot_dir.'/PNG';
		my $plot_dir_svg = $plot_dir.'/SVG';

		for my $dir ($plot_dir,$plot_dir_png,$plot_dir_svg){
			unless (-d $dir){
				mkdir($dir,0755) or die "Can't create $dir: !$\n";
			}
		}

		for my $orientation ('normal', 'inverted'){

			my $conf = $circos_cat_dir."/concatenated.$affix.$orientation.conf";
			my $circos_plot = "$circos_prefix.$affix.$orientation.png";

			$circos_todo_list{$conf}{'dir'} = $plot_dir;
			$circos_todo_list{$conf}{'png'} = $circos_plot;

		}

	}

	elsif ($mode eq 'pair'){

		my $plot_dir = $circos_plot_dir.'/Pairwise';
		my $plot_dir_png = $plot_dir.'/PNG';
		my $plot_dir_svg = $plot_dir.'/SVG';

		for my $dir ($plot_dir,$plot_dir_png,$plot_dir_svg){
			unless (-d $dir){
				mkdir($dir,0755) or die "Can't create $dir: !$\n";
			}
		}

		opendir(PAIRS, $circos_pair_dir) or die "Can't open $circos_pair_dir: $!\n";
		my @pairs;
		while (my $entry = readdir(PAIRS)){
			if ($entry =~ /_vs_/){
				push(@pairs, $entry);
			}
		}

		for my $pair (@pairs){

			for my $orientation ('normal', 'inverted'){

				my $conf = $circos_pair_dir.'/'.$pair.'/'.$pair.'.'.$affix.'.'.$orientation.'.conf';
				my $circos_plot = "$pair.$affix.$orientation.png";

				$circos_todo_list{$conf}{'dir'} = $plot_dir;
				$circos_todo_list{$conf}{'png'} = $circos_plot;

			}

		}

	}

}

sub run_circos {

	if ($circos_orientation ne 'both'){
		$circos_num = $circos_num / 2;
	}

	while (my $conf = shift@circos_files){

		if ($circos_orientation eq 'normal'){
			if ($conf =~ /\.inverted\.conf$/){
				next;
			}
		}
		elsif ($circos_orientation eq 'inverted'){
			if ($conf =~ /\.normal\.conf$/){
				next;
			}
		}

		my $plot_dir = $circos_todo_list{$conf}{'dir'};
		my $circos_plot = $circos_todo_list{$conf}{'png'};

		## Progess counter
		my $x;
		if ($circos_orientation eq 'normal'){
			$x = $circos_num - (scalar(@circos_files) / 2);
		}
		elsif ($circos_orientation eq 'inverted'){
			$x = $circos_num - ((scalar(@circos_files) + 1) / 2);
		}
		else {
			$x = $circos_num - scalar(@circos_files);

		}
		print "$x / $circos_num - plotting $plot_dir/$circos_plot\n";

		## Plot circos configuration file
		system ("
			circos \\
				-conf $conf \\
				-outputdir $plot_dir \\
				-outputfile $circos_plot \\
				-silent
		");

	}
	threads->exit();
}
