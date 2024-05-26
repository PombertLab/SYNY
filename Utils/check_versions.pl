#!/usr/bin/env perl
## Pombert Lab, 2024

my $name = 'check_versions.pl';
my $version = '0.1b';
my $updated = '2020-05-22';

use strict;
use warnings;
use File::Find;
use File::Basename;
use Cwd qw(abs_path);
use Getopt::Long qw(GetOptions);

my $usage = <<"USAGE";
NAME        ${name}
VERSION     ${version}
UPDATED     ${updated}
SYNOPSIS    Checks SYNY and subscripts versions

COMMAND     ${name} \\
              --dir ./ \\
              --out versions.txt

OPTIONS:
-d (--dir)      Directory to check
-o (--out)      Output text file [Default: versions.txt]
-g (--git)      Specify git tag (otherwise finds latest from .git/FETCH_HEAD)
USAGE

unless (@ARGV){
	print "\n$usage\n";
	exit(0);
};

my $indir;
my $outfile = 'versions.txt';
my $git_tag;
GetOptions(
    'd|dir=s' => \$indir,
    'o|out=s' => \$outfile,
    'g|git=s' => \$git_tag
);

## Locating scripts/subscripts
my @scripts;
$indir = abs_path($indir);
find(\&scripts, $indir);

# Getting git tag automatically
unless ($git_tag){
    my $gitfile = $indir.'/.git/FETCH_HEAD';
    open GIT, '<', $gitfile or die "Can't read $gitfile: $!\n";
    while (my $line = <GIT>){
        chomp $line;
        my @data = split("\t", $line);
        my $tag_id = $data[2];
        if ($tag_id =~ /tag \'(\S+)\'/){
            $git_tag = $1;
        }
    }
}

## Getting scripts + versions info
my %master_scripts;
my %subscripts;

for my $abscript (@scripts){

    my ($script,$path) = fileparse($abscript);
    my $subdir;
    my $script_version;

    ## Script location
    if ($path =~ /\/(\w+)\/$/){
        my $tmp = $1;
        if ($tmp =~ /SYNY/){
            $subdir = 'Main';
        }
        else {
            $subdir = $tmp;
        }
    }

    ## Grab version from script
    open SC, "<", $abscript or die "Can't read $abscript: $!\n";
    while (my $line = <SC>){

        chomp $line;

        if ($line =~ /^(my \$)?version = \'(\S+)\'/){

            $script_version = $2;

            if ($subdir eq 'Main'){
                $master_scripts{$script} = $script_version;
            }
            else{
                $subscripts{$subdir}{$script} = $script_version;
            }

        }
    }

}

## SYNY info
open OUT, '>', $outfile or die "Can't create $outfile: $!\n";
for my $fh (\*OUT, \*STDOUT){

    print $fh "\n";
    print $fh '#############################################################'."\n";
    print $fh 'SYNY - Version '.$git_tag."\n";
    print $fh '#############################################################'."\n";

    ## Main dir
    my $header = '# Main';
    my $hl = length($header);
    my $pad = 40 - $hl;

    print $fh "\n".$header.' ' x $pad.'Versions'."\n";

    for my $subscript (sort(keys %master_scripts)){
        $hl = length($subscript);
        $pad = 40 - $hl;
        print $fh $subscript.' ' x $pad.$master_scripts{$subscript}."\n";
    }

    ## Subdirs
    for my $subdir (sort(keys %subscripts)){

        print $fh "\n".'# '.$subdir."\n";

        for my $subscript (sort(keys %{$subscripts{$subdir}})){
            $hl = length($subscript);
            $pad = 40 - $hl;
            print $fh $subscript.' ' x $pad.$subscripts{$subdir}{$subscript}."\n";
        }

    }

    print $fh "\n";
    print $fh '#############################################################'."\n";

}

### Subroutine(s)
sub scripts {

    my $sub = abs_path($_);

    if (($sub =~ /.pl$/) or ($sub =~ /.py$/)){
        push (@scripts, $sub);
    }

}