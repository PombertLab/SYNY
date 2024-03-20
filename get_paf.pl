#!/usr/bin/env perl
## Pombert Lab, 2024

use strict;
use warnings;
use File::Basename;
use Cwd qw(abs_path);
use Getopt::Long qw(GetOptions);

my $name = 'get_paf.pl';
my $version = '0.3c';
my $updated = '2024-03-20';

my $usage =<<"USAGE";
NAME        ${name}
VERSION     ${version}
UPDATED     ${updated}
SYNOPSIS    Performs pairwize genome colinearity (paf,maf,aln) comparisons using Minimap2

EXAMPLE     ${name} \\
              -f *.fasta \\
              -o PAFR_5 \\
              -resume \\
              -asm 5

PREREQS     Minimap2:     https://github.com/lh3/minimap2

OPTIONS:
-f (--fasta)    FASTA files to compare
-o (--outdir)   Output directory [Default: ./PAF]
-r (--resume)   Resume computation (skip completed alignments)
-a (--asm)      Specify minimap2 max divergence preset (asm 5, 10 or 20) [Default: off]
USAGE

die "\n$usage\n" unless @ARGV;
my @commands = @ARGV;

my @fasta;
my $outdir = './PAF';
my $resume;
my $asm;
GetOptions(
    'f|fasta=s@{1,}' => \@fasta,
    'o|outdir=s' => \$outdir,
    'r|resume' => \$resume,
    'asm=i' => \$asm
);

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

push (@subdirs, ($paf_dir,$maf_dir,$blast_dir));

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

foreach my $query (@fasta){

    my $bquery = fileparse($query);
    $bquery =~ s/\.(\w+)$//;

    foreach my $target (@fasta){

        my $btarget = fileparse($target);
        $btarget =~ s/\.(\w+)$//;

        unless ($query eq $target){

            print STDOUT "\n".'Minimap2: '.$bquery.' vs. '.$btarget."\n\n";
            print LOG 'Minimap2: '.$bquery.' vs. '.$btarget;

            # Minimap2 PAF, MAF and ALN (BLAST-like) output files
            my $tmp_paf_outfile = $paf_dir.'/'.$bquery.'_vs_'.$btarget.'.tmp.paf';
            my $paf_outfile = $paf_dir.'/'.$bquery.'_vs_'.$btarget.'.paf';
            my $maf_outfile = $maf_dir.'/'.$bquery.'_vs_'.$btarget.'.maf';
            my $blast_outfile = $blast_dir.'/'.$bquery.'_vs_'.$btarget.'.aln';

            my $map_time_start = time;

            # Skip alignment if paf file is found
            if ((-e $paf_outfile) && ($resume)){
                print "Found $paf_outfile, skipping minimap2 alignment ...\n";
                goto GETMAF;
            }

            # Running minimap2 (PAF output)
            system(
                "minimap2 \\
                    $asm_flag \\
                    -c \\
                    --cs=long \\
                    $query \\
                    $target \\
                    > $tmp_paf_outfile
                "
            ) == 0 or checksig();

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