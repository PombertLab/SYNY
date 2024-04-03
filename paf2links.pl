#!/usr/bin/env perl
## Pombert Lab, 2024

my $name = 'paf2links.pl';
my $version = '0.1';
my $updated = '2024-04-03';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;
use File::Path qw(make_path);

my $usage = <<"USAGE";
NAME        ${name}
VERSION     ${version}
UPDATED     ${updated}
SYNOPSIS    Creates links for Circos from PAF files

COMMAND     ${name} \\
              -paf PAF_DIR \\
              -links paf_links.txt

OPTIONS:
-p (--paf)          PAF file(s) directory
-l (--links)        Desired Circos output link file [Default: paf_links.txt]
--clusters          Color by cluster [Default: off]
--custom_file       Load custom colors from file
--custom_preset     Use a custom color preset; e.g.
                    # chloropicon - 20 colors - Lemieux et al. (2019) https://pubmed.ncbi.nlm.nih.gov/31492891/
                    # encephalitozoon - 11 colors - Pombert et al. (2012) https://pubmed.ncbi.nlm.nih.gov/22802648/
USAGE
die "\n$usage\n" unless @ARGV;

my $paf_dir;
my $links_file = 'paf_links.txt';
my $clusters;
my $custom_file;
my $custom_cc;
GetOptions(
    'p|paf=s' => \$paf_dir,
    'l|links=s' => \$links_file,
    'clusters' => \$clusters,
    'custom_file=s' => \$custom_file,
    'custom_preset=s' => \$custom_cc
);

## Grabbing paf files
opendir (PAFDIR, $paf_dir) or die "\n\n[ERROR]\tCan't open $paf_dir: $!\n\n";

my @paf_files;
while (my $file = readdir(PAFDIR)){
	if ($file =~ /\.paf$/){
		push (@paf_files, "$paf_dir/$file");
	}
}

#################### Circos colors
my @yellows = ('vlyellow','lyellow','yellow','dyellow','vdyellow','vvdyellow');
my @reds = ('vvlred','vlred','lred','red','dred','vdred','vvdred');
my @oranges = ('vvlorange','vlorange','lorange','orange','dorange','vdorange','vvdorange');
my @greens = ('vvlgreen','vlgreen','lgreen','green','dgreen','vdgreen','vvdgreen');
my @blues = ('vvlblue','vlblue','lblue','blue','dblue','vdblue','vvdblue');
my @purples = ('vvlpurple','vlpurple','lpurple','purple','dpurple','vdpurple','vvdpurple');
my @greys = ('vvlgrey','vlgrey','lgrey','grey','dgrey','vdgrey','vvdgrey','vvvdgrey');

my @rainbow = (@reds,@oranges,@yellows,@greens,@blues,@purples); ## 41 colors total
my @bowgrey = (@rainbow,@greys); ## 49 colors total

my @color_set = @bowgrey;

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

## Creating Circos links file
open PLINK, '>', $links_file or die "Can't create $links_file: $!\n";
print PLINK '#locus1 start end locus2 start end'."\n";

my $zdepth = 1;
my $link_color = 'grey_a5';
my $color_number = 0;

for my $paf_file (@paf_files){

	open PAF, '<', $paf_file or die "Can't read $paf_file: $!\n";

	$color_number = 0;

	while (my $line = <PAF>){

		chomp $line;

		if ($clusters){
			$link_color = $color_set[$color_number];
			$color_number++;
			if ($color_number >= scalar(@color_set)){
				$color_number = 0;
			}
		}

		my @data = split("\t", $line);
		my $locus1 = $data[0];
		my $l1_start = $data[2];
		my $l1_end = $data[3];
		my $locus2 = $data[5];
		my $l2_start = $data[7];
		my $l2_end = $data[8];

		print PLINK "$locus1 $l1_start $l1_end $locus2 $l2_start $l2_end ";
		print PLINK 'color='.$link_color;
		print PLINK ',z='.$zdepth."\n";

	}
	close PAF;
}
close PLINK;

####################################################################
### Subroutine(s)
####################################################################

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