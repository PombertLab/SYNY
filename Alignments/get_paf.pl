#!/usr/bin/env perl
## Pombert Lab, 2024

my $name = 'get_paf.pl';
my $version = '0.4b';
my $updated = '2024-05-27';

use strict;
use warnings;
use File::Basename;
use Cwd qw(abs_path);
use Getopt::Long qw(GetOptions);

#########################################################################
### Command line options
#########################################################################

my $usage =<<"USAGE";
NAME        ${name}
VERSION     ${version}
UPDATED     ${updated}
SYNOPSIS    Performs pairwize genome colinearity (paf,maf,aln) comparisons using Minimap2

EXAMPLE     ${name} \\
              -f *.fasta \\
              -o PAFR_5 \\
              -a minimap \\
              -threads 8 \\
              -resume \\
              -asm 5

PREREQS     Minimap2:     https://github.com/lh3/minimap2

OPTIONS:
-f (--fasta)    FASTA files to compare
-o (--outdir)   Output directory [Default: ./PAF]
-a (--aligner)  Genome alignment tool, minimap2 or mashmap3 [Default: minimap2]
-r (--resume)   Resume computation (skip completed alignments)
-t (--threads)  Number of threads for minimap2 [Default: 8]
-a (--asm)      Specify minimap2 max divergence preset (asm 5, 10 or 20) [Default: off]
-p (--percent)  Specify mashmap3 percentage identity [Default: 75] 
-v (--version)  Show script version
USAGE

unless (@ARGV){
	print "\n$usage\n";
	exit(0);
};

my @commands = @ARGV;

my @fasta;
my $outdir = './PAF';
my $aligner = 'minimap2';
my $resume;
my $threads = 8;
my $asm;
my $mashmap_pid = 75;
my $sc_version;
GetOptions(
    'f|fasta=s@{1,}' => \@fasta,
    'o|outdir=s' => \$outdir,
    'a|aligner=s' => \$aligner,
    't|threads=i' => \$threads,
    'r|resume' => \$resume,
    'asm=i' => \$asm,
    'p|percent=s' => \$mashmap_pid,
    'v|version' => \$sc_version
);

#########################################################################
### Version
#########################################################################

if ($sc_version){
    print "\n";
    print "Script:     $name\n";
    print "Version:    $version\n";
    print "Updated:    $updated\n\n";
    exit(0);
}

#########################################################################
### Output dir/subdirs
#########################################################################

# Main dir
$outdir = abs_path($outdir);
my @subdirs = ($outdir);

# Subdirs
my $paf_dir = $outdir.'/'.'PAF';
my $maf_dir = $outdir.'/'.'MAF';
my $blast_dir = $outdir.'/'.'ALN';

push (@subdirs, $paf_dir);

if ($aligner =~ /minimap/){
    push (@subdirs, ($maf_dir,$blast_dir));
}

for my $dir (@subdirs){
    unless (-d $dir){
        mkdir ($dir, 0755) or die "Can't create $dir: $!";
    }
}

#########################################################################
### Creating log
#########################################################################

my $start = localtime();
my $tstart = time;
my $date = `date`;

my $logfile = $outdir.'/paf.log';
open LOG, '>', $logfile or die "Can't create file $logfile: $!\n";

print LOG "\# $date\n";
print LOG "COMMAND LINE:\n";
print LOG "$0 @commands\n\n";

#########################################################################
### Align genomes with minimap2
#########################################################################

my $asm_flag = '';
if ($asm){
    $asm_flag = '-x asm'.$asm;
}

my $fasta_num = scalar(@fasta);
my $to_do = $fasta_num * ($fasta_num - 1);
my $current_iteration = 0;

foreach my $query (@fasta){

    my $bquery = fileparse($query);
    $bquery =~ s/\.(\w+)$//;

    foreach my $target (@fasta){

        my $btarget = fileparse($target);
        $btarget =~ s/\.(\w+)$//;

        unless ($query eq $target){

            # $current_iteration++;

            print STDOUT "\n"."$aligner: ".++$current_iteration.'/'.$to_do.' --- '.$bquery.' vs. '.$btarget."\n\n";
            print LOG "$aligner: ".$bquery.' vs. '.$btarget;

            my $affix = 'mmap';

            # Minimap2 PAF, MAF and ALN (BLAST-like) output files
            my $tmp_paf_outfile = $paf_dir.'/'.$bquery.'_vs_'.$btarget.'.tmp.paf';
            my $paf_outfile = $paf_dir.'/'.$bquery.'_vs_'.$btarget.".$affix.paf";
            my $maf_outfile = $maf_dir.'/'.$bquery.'_vs_'.$btarget.".$affix.maf";
            my $blast_outfile = $blast_dir.'/'.$bquery.'_vs_'.$btarget.".$affix.aln";

            my $map_time_start = time;

            # Skip alignment if paf file is found
            if ((-e $paf_outfile) && ($resume)){
                print "Found $paf_outfile, skipping $aligner alignment ...\n";
                goto GETMAF;
            }

            if ($aligner =~ /minimap/i){
                # Running minimap2 (PAF output)
                system(
                    "minimap2 \\
                        -t $threads \\
                        $asm_flag \\
                        -c \\
                        --cs=long \\
                        $query \\
                        $target \\
                        > $tmp_paf_outfile
                    "
                ) == 0 or checksig();
            }

            elsif ($aligner =~ /mashmap/i){
                # Running mashmap3 (PAF output)
                system(
                    "mashmap \\
                        -t $threads \\
                        -q $target \\
                        -r $query \\
                        --perc_identity $mashmap_pid \\
                        -o $tmp_paf_outfile
                    "
                ) == 0 or checksig();
            }

            ## Sorting PAF output by genomic coordinates
            ## PAF files can be large, so keeping only the sorted version
            ## to save on disk usage

            open PAF, '<', $tmp_paf_outfile or die "Can't open $tmp_paf_outfile: $!";
            open SPAF, '>', $paf_outfile or die "Can't create $paf_outfile: $!\n";

            my %db;

            while (my $line = <PAF>){

                chomp $line;
                my @data = split("\t", $line);
                my $locus = $data[0];
                my $start = $data[2];

                ## Pushing as an array in case multiple matches start at the same position
                push (@{$db{$locus}{$start}}, $line);

            }

            foreach my $locus (sort (keys %db)){
                foreach my $start (sort {$a <=> $b}(keys %{$db{$locus}}) ){
                    my @tmp = @{$db{$locus}{$start}};
                    foreach my $line (@tmp){
                        print SPAF "$line\n";
                    }
                }
            }

            close PAF;
            close SPAF;

            system ("rm $tmp_paf_outfile");

            ## Creating MAF file
            GETMAF:

            if ($aligner =~ /minimap/i){

                if ((-e $maf_outfile) && ($resume)){
                    print "Found $maf_outfile, skipping paftools.js conversion ...\n";
                    goto GETALN;
                }

                system (
                    "paftools.js view \\
                    -f maf \\
                    $paf_outfile \\
                    > $maf_outfile
                    "
                ) == 0 or checksig();

                ## Creating BLAST file
                GETALN:
                if ((-e $blast_outfile) && ($resume)){
                    print "Found $blast_outfile, skipping paftools.js conversion ...\n";
                    goto ENDLOG;
                }

                system (
                    "paftools.js view \\
                    -f aln \\
                    $paf_outfile \\
                    > $blast_outfile
                    "
                ) == 0 or checksig();
            }

            ENDLOG:
            my $map_time_end = time;
            my $map_time_taken = $map_time_end - $map_time_start;
            print LOG "; $map_time_taken seconds\n";

        }
    }
}

###################################################################################################
## Subroutines
###################################################################################################

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