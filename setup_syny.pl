#!/usr/bin/env perl
# Pombert Lab, 2024

my $name = 'setup_syny.pl';
my $version = '0.3';
my $updated = '2024-05-07';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;
use Cwd qw(abs_path);
use Cwd qw(getcwd);
use File::Path qw(make_path);

my $usage = <<"USAGE";
NAME        ${name}
VERSION     ${version}
UPDATED     ${updated}
SYNOPSIS    Install SYNY dependencies and adds install locations to the 
            specified configuration file

EXAMPLE     setup_SYNY.pl \\
              -l fedora \\
              -c test.cfg \\
              -i ./TOOLS

OPTIONS:
-l (--linux)    Linux distribution [Default: fedora]
-c (--config)   Configuration file to edit/create [Default: ~/.bash_profile]
-i (--install)  Installation directory (for dependencies) [Default: ./TOOLS]
--list_distro   Lists supported distribution/package managers
USAGE
die "\n$usage\n" unless @ARGV;

my $linux = 'fedora';
my $config = '~/.bash_profile';
my $install_dir = './TOOLS';
my $list_distro;
GetOptions(
    'l|linux=s' => \$linux,
    'c|config=s' => \$config,
    'i|install=s' => \$install_dir,
    'list_distro' => \$list_distro
);

###################################################################################################
## Grabbing absolute paths
###################################################################################################

my $init_dir = getcwd();
$install_dir = abs_path($install_dir);
$config = abs_path($config);
my ($script,$syny_path) = fileparse(abs_path($0));

###################################################################################################
## Installing Linux dependencies (requires sudo)
###################################################################################################

### Linux distro/package manager check
my %linux_distros = (
    'fedora' => 'dnf',
    'debian' => 'apt',
    'ubuntu' => 'apt',
    'kali' => 'apt',
    'opensuse' => 'zypper',
);

$linux = lc($linux);

sub distros{

    print "\nSupported Linux distributions/package managers are:\n\n";

    foreach my $distro (sort(keys %linux_distros)){
        my $slen = 12 - length($distro);
        my $spacer = ' ' x $slen;
        print $distro.$spacer.$linux_distros{$distro}."\n";
    }

    print "\n";

}

unless (exists $linux_distros{$linux}){

    print "\nUnrecognized Linux distribution: $linux\n";
    distros();
    print "Exiting...\n\n";
    exit();

}

if ($list_distro){
    distros();
    exit();
}

### Installing deps.

print "\nPlease enter sudo password to install $linux dependencies\n\n";

my $package_manager = $linux_distros{$linux};

if ($package_manager eq 'apt'){
    system ("
      sudo \\
        apt install -y \\
        git \\
        curl \\
        build-essential \\
        zlib1g-dev \\
        python3-matplotlib \\
        python3-seaborn \\
        python3-pandas \\
        python3-scipy \\
        libperlio-gzip-perl \\
        cpanminus \\
        libgd-perl
    ") == 0 or checksig();
}

elsif ($package_manager eq 'dnf'){

    system ('sudo dnf group install -y "C Development Tools and Libraries"');

    system ("
        sudo \\
        dnf install -y \\
        git \\
        curl \\
        zlib-devel \\
        python3-matplotlib \\
        python3-seaborn \\
        python3-pandas \\
        python3-scipy \\
        perl-PerlIO-gzip \\
        perl-App-cpanminus \\
        perl-GD
    ") == 0 or checksig();

}

elsif ($package_manager eq 'zypper'){

    system ("
        sudo \\
        zypper install -y \\
        git \\
        curl \\
        patterns-devel-base-devel_basis \\
        zlib-devel \\
        python3-matplotlib \\
        python3-seaborn \\
        python3-pandas \\
        python3-scipy \\
        perl-PerlIO-gzip \\
        perl-App-cpanminus \\
        perl-GD
    ") == 0 or checksig();

}

print "\nInstalling Circos perl dependencies...\n\n";

system ("
    sudo cpanm \\
        Roman \\
        Clone \\
        Config::General \\
        Font::TTF::Font \\
        Getopt::ArgvFile \\
        List::MoreUtils \\
        Math::Bezier \\
        Math::Round \\
        Math::VecStat \\
        Params::Validate \\
        Readonly \\
        Regexp::Common \\
        Set::IntSpan \\
        Statistics::Basic \\
        SVG \\
        Text::Format
") == 0 or checksig();

###################################################################################################
## Creating installation directory
###################################################################################################

unless(-d $install_dir){
	make_path($install_dir,{mode => 0755}) or die "Can't create $install_dir: $!\n";
}

###################################################################################################
## Writing/appending to configuration file
###################################################################################################

my $diamond = '>';
if (-e $config){
    $diamond = '>>';
}

open CFG, "$diamond", $config or die "Can't open $config: $!\n";
print CFG "\n\n";
print CFG 'PATH=$PATH:'.$syny_path.'                    ## SYNY'."\n";
print CFG 'PATH=$PATH:'.$syny_path.'Alignments          ## SYNY'."\n";
print CFG 'PATH=$PATH:'.$syny_path.'Clusters            ## SYNY'."\n";
print CFG 'PATH=$PATH:'.$syny_path.'Examples            ## SYNY'."\n";
print CFG 'PATH=$PATH:'.$syny_path.'Plots               ## SYNY'."\n";
print CFG 'PATH=$PATH:'.$syny_path.'Utils               ## SYNY'."\n";


###################################################################################################
## Installing Circos
###################################################################################################

print "\nInstalling Circos...\n\n";

chdir $install_dir;

system("
    curl \\
    -L https://circos.ca/distribution/circos-0.69-9.tgz \\
    --insecure \\
    -o $install_dir/circos-0.69-9.tgz
") == 0 or checksig();

system("tar -zxvf $install_dir/circos-0.69-9.tgz --directory $install_dir");
system("rm $install_dir/circos-0.69-9.tgz");

my $circos_dir = $install_dir.'/circos-0.69-9/bin';
print CFG 'PATH=$PATH:'.$circos_dir.'         ## Circos'."\n";

###################################################################################################
## Installing Diamond
###################################################################################################

print "\nInstalling Diamond...\n\n";

chdir $install_dir;

my $diamond_dir = $install_dir.'/diamond';
my $diamond_ver = 'v2.1.9';

unless(-d $diamond_dir){
	make_path($diamond_dir,{mode => 0755}) or die "Can't create $diamond_dir: $!\n";
}

system ("
    curl \\
    -L https://github.com/bbuchfink/diamond/releases/download/$diamond_ver/diamond-linux64.tar.gz \\
    -o $install_dir/diamond-linux64.tar.gz
") == 0 or checksig();

system("
    tar -zxvf $install_dir/diamond-linux64.tar.gz --directory $diamond_dir
    rm $install_dir/diamond-linux64.tar.gz
");

print CFG 'PATH=$PATH:'.$diamond_dir.'                   ## Diamond'."\n";

###################################################################################################
## Installing Minimap2
###################################################################################################

print "\nInstalling Minimap2...\n\n";

chdir $install_dir;
my $minimap_dir = $install_dir.'/minimap2';

if (-d $minimap_dir){
    system("rm -R -f $minimap_dir");
}

system("git clone https://github.com/lh3/minimap2 $minimap_dir");
chdir $minimap_dir;
system("make") == 0 or checksig();

system("
    curl \\
    -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 \\
    -o $install_dir/k8-0.2.4.tar.bz2
") == 0 or checksig();

system("tar -jxf $install_dir/k8-0.2.4.tar.bz2 --directory $minimap_dir");
system("cp $minimap_dir/k8-0.2.4/k8-`uname -s` $minimap_dir/k8");
system("rm -R $minimap_dir/k8-0.2.4/");
system("rm $install_dir/k8-0.2.4.tar.bz2");

print CFG 'PATH=$PATH:'.$minimap_dir.'                  ## Minimap2'."\n";
print CFG 'PATH=$PATH:'.$minimap_dir.'/misc             ## Minimap2 paftools.js'."\n";

###################################################################################################
## Completion
###################################################################################################

print CFG "\n";
print CFG 'export PATH'."\n";
close CFG;

print "\nCompleted. Please run source on your configuration file, i.e.:\n\n";
print "source $config\n\n";

###################################################################################################
## Subroutine(s)
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