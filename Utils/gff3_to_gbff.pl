#!/usr/bin/env perl
## Pombert Lab, 2024

my $name = 'gff3_to_gbff.pl';
my $version = '0.2b';
my $updated = '2026-02-04';

use strict;
use warnings;
use File::Basename;
use File::Path qw(make_path);
use PerlIO::gzip;
use Getopt::Long qw(GetOptions);

#########################################################################
### Command line options
#########################################################################

my $usage =<<"USAGE";
NAME        ${name}
VERSION     ${version}
UPDATED     ${updated}
SYNOPSIS    Converts GFF3 + FASTA files to pseudo-GenBank flat file (GBFF) format
            for use with SYNY

REQS        GFF3 and FASTA files must have the same file prefixes; e.g.
            genome_1.fna genome_1.gff

NOTE        Tested with NCBI/AGAT GFF3/GTF files;
            expects gene/(mRNA|transcript)/exon/CDS entries in the 'type' column

COMMAND     $name \\
              --fasta *.fasta \\
              --gff3 *.gff \\
              --outdir ./ \\
              --gcode 1 \\
              --gzip

OPTIONS:
-f (--fasta)    FASTA file(s) to convert (gziped files are supported)
-g (--gff3)     GFF3 files to convert (gziped files are supported)
-t (--type)     GFF3 type (ncbi, agat) [Default: ncbi]
-o (--outdir)   Output directory [Default: GBFF]
-z (--gzip)     Compress the GBFF output files
-i (--id)       Use the ID field as /product ## Can be useful if GFF3 lacks product entries
-c (--gcode)    NCBI genetic code [Default: 1]
                1  - The Standard Code
                2  - The Vertebrate Mitochondrial Code
                3  - The Yeast Mitochondrial Code
                4  - The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
                11 - The Bacterial, Archaeal and Plant Plastid Code
                NOTE - For complete list; see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
-v (--verbose)  Add verbosity
--version       Show script version
USAGE

unless (@ARGV){
    print "\n$usage\n";
    exit(0);
};

my @fasta;
my @gff3;
my $gtype = 'ncbi';
my $outdir = 'GBFF';
my $gzip_flag;
my $id_flag;
my $gc = 1;
my $verbose;
my $sc_version;
GetOptions(
    'f|fasta=s@{1,}' => \@fasta,
    'g|gff3=s@{1,}' => \@gff3,
    't|type=s' => \$gtype,
    'o|outdir=s' => \$outdir,
    'z|gzip' => \$gzip_flag,
    'i|id' => \$id_flag,
    'c|gcode=i' => \$gc,
    'v|verbose' => \$verbose,
    'version' => \$sc_version
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
### Checking gff3 type
#########################################################################

my %known_types = (
    'ncbi' => '1',
    'agat' => '1'
);

my $entry = lc($gtype);
unless (exists $known_types{$entry}){
    print ("\nFile type $entry is unknown. Possible entries are: ");
    foreach my $key (sort (keys %known_types)){
        print $key." ";
    }
    print "\n\n";
    exit(0);
}

#########################################################################
### Pigz check
#########################################################################

my $gzip_tool = 'gzip';
my $pigz_check = `echo \$(command -v pigz)`;
if ($pigz_check =~ /pigz/){
    $gzip_tool = 'pigz';
}

#########################################################################
### Output dir/subdirs
#########################################################################

unless (-d $outdir){
    make_path($outdir,{mode=>0755}) or die "Can't create $outdir: $!\n";
}

#########################################################################
### Loading NCBI genetic codes
#########################################################################
my %gcodes;
gcodes();

#########################################################################
### Grabbing basename from GFF3 files
#########################################################################

my %gff3;

while (my $gff3 = shift@gff3){
    my ($basename,$path) = fileparse($gff3);
    $basename =~ s/\.gz$//;
    $basename =~ s/\.\w+$//;
    $gff3{$basename} = $gff3;
}

#########################################################################
### Parsing fasta/gff3 file(s)
#########################################################################

while (my $fasta = shift@fasta){

    if ($verbose){
        print 'Parsing           : '.$fasta."\n";
    }

    ## Checking for gzip file extension
    my $gzip = '';
    if ($fasta =~ /.gz$/){
        $gzip = ':gzip';
    }

    ## Grabbing basename + I/O
    open FASTA, "<$gzip", $fasta or die "Can't read $fasta: $!\n";
    my ($basename,$path) = fileparse($fasta);
    $basename =~ s/\.gz$//;
    $basename =~ s/\.\w+$//;
    my $outfile = $outdir.'/'.$basename.'.gbff';
    open GBFF, '>', $outfile or die "Can't create $outfile: $!\n";

    ## Creating a database of sequences
    my %sequences;
    my %contigs;
    my $header;

    while (my $line = <FASTA>){

        chomp $line;
        if ($line =~ /^>(\w+)/){
            $header = $1;
            @{$contigs{$header}} = ();
        }
        else {
            $sequences{$header} .= $line;
        }

    }

    ## Linking gff3 to fasta file using its file prefix (must be the same)
    my $gff3 = $gff3{$basename};
    my $gz = '';
    if ($gff3 =~ /.gz$/){
        $gz = ':gzip';
    }
    open GFF3, "<$gz", $gff3 or die "Can't read $gff3: $!\n";
    
    if ($verbose){
        print 'Parsing           : '.$gff3."\n";
    }

    ## Parsing GFF3 file
    my %genes;
    my %rnas;
    my $gene_counter = 0;

    my %agat_products;

    while (my $line = <GFF3>){

        chomp $line;

        # Skip comments
        if ($line =~ /^#/){
            next;
        }

        # Skip blank lines
        elsif ($line =~ /^\s*$/){
            next;
        }

        else{

            my @data = split("\t", $line);

            my ($seqid) = $data[0] =~ /^(\w+)/;
            my $source = $data[1];
            my $type = $data[2];
            my $start = $data[3];
            my $end = $data[4];
            my $score = $data[5];
            my $strand = $data[6];
            my $phase = $data[7];
            my @attributes = split(";", $data[8]);

            my $feature_id;
            my $parent_id;
            my $feature_name;
            my $locus_tag;
            my $product;

            foreach my $attribute (@attributes){

                if ($entry eq 'ncbi'){
                    if ($attribute =~ /^ID=(.*)/){
                        $feature_id = $1; 
                    }
                    if ($attribute =~ /^Parent=(.*)/){
                        $parent_id = $1; 
                    }
                    if ($attribute =~ /^Name=(.*)/){
                        $feature_name = $1; 
                    }
                    if ($attribute =~ /^locus_tag=(.*)/){
                        $locus_tag = $1; 
                    }
                    if ($attribute =~ /^product=(.*)/){
                        $product = $1; 
                    }
                }

                elsif ($entry eq 'agat'){
                    if ($attribute =~ /^\s?gene_id \"(.*)\"/){
                        $locus_tag = $1;
                        if ($type eq 'gene'){
                            $feature_name = $locus_tag;
                        }
                    }
                    if ($attribute =~ /^\s?ID \"(.*)\"/){
                        $feature_id = $1;
                        if  ($type ne 'gene'){
                            $feature_name = $feature_id;
                        }
                    }
                    if ($attribute =~ /^\s?Name \"(.*)\"/){
                        $product = $1;
                        $agat_products{$feature_id} = $product;
                    }
                    if ($attribute =~ /^\s?Parent \"(.*)\"/){
                        $parent_id = $1; 
                    }
                }

            }

            ## Renaming agat transcript key to mRNA
            if ($entry eq 'agat'){
                if ($type eq 'transcript'){
                    $type = 'mRNA';
                }
            }

            if ($type eq 'gene'){
                
                ## Keeping track of the features per contig
                push (@{$contigs{$seqid}}, $feature_id);
                
                $gene_counter++;
                unless ($locus_tag){
                    $locus_tag = $seqid.'_'.$gene_counter;
                }
                $genes{$feature_id}{'locus'} = $locus_tag;
                $genes{$feature_id}{'strand'} = $strand;
                $genes{$feature_id}{'start'} = $start;
                $genes{$feature_id}{'end'} = $end;
            
            }

            elsif ($type =~ /RNA/){
                push(@{$genes{$parent_id}{'children'}}, $feature_id);
                $genes{$parent_id}{'type'} = $type;
                $rnas{$feature_id}{'type'} = $type;
                $rnas{$feature_id}{'strand'} = $strand;
                $rnas{$feature_id}{'parent'} = $parent_id;
            }

            elsif ($type eq 'exon'){

                my $coordinates = $start.'..'.$end;
                push(@{$rnas{$parent_id}{'coordinates'}}, $coordinates);

                ## AGAT doesn't put products on exon lines (in the files tested)
                if ($entry eq 'agat'){
                    if (exists $agat_products{$parent_id}){
                        $rnas{$parent_id}{'product'} = $agat_products{$parent_id};
                    }
                    else{
                        $rnas{$parent_id}{'product'} = 'undefined product';
                    }
                }
                else{
                    if ($product){
                        $rnas{$parent_id}{'product'} = $product;
                    }
                    else{
                        $rnas{$parent_id}{'product'} = 'undefined product';
                    }
                }
            }

            elsif ($type eq 'CDS'){

                my $coordinates = $start.'..'.$end;
                push(@{$rnas{$parent_id}{'coordinates_cds'}}, $coordinates);

                ## AGAT doesn't put products on CDS lines (in the files tested)
                if ($entry eq 'agat'){
                    if (exists $agat_products{$parent_id}){
                        $rnas{$parent_id}{'product'} = $agat_products{$parent_id};
                    }
                    else{
                        $rnas{$parent_id}{'product'} = 'undefined product';
                    }
                }

                else{
                    if ($product){
                        $rnas{$parent_id}{'product'} = $product;
                    }
                    else{
                        $rnas{$parent_id}{'product'} = 'undefined product';
                    }
                }

            }

        }
    }

    ## Iterating though each sequence found in FASTA file
    for my $sequence (sort(keys %sequences)){

        ## GBFF header
        my $length = length $sequences{$sequence};
        print GBFF "LOCUS       $sequence           $length bp    DNA"."\n";
        print GBFF "DEFINITION  $sequence"."\n";
        print GBFF "ACCESSION   $sequence"."\n";
        print GBFF "VERSION     $sequence"."\n";
        print GBFF "COMMENT     Pseudo-GBFF file for SYNY"."\n"; 
        print GBFF "KEYWORDS    ."."\n";
        print GBFF "SOURCE      ."."\n";
        print GBFF "ORGANISM    ."."\n";
        print GBFF "            ."."\n";
        print GBFF "FEATURES             Location/Qualifiers"."\n";

        ## Annotations
        my @features = @{$contigs{$sequence}};

        ## Sorting annotations by start positions => annotations in GFF3 files
        ## can be out-of-order

        my %feat;
        foreach my $feature (@features){
            if (exists $feat{$feature}){
                print "ID = $feature is not unique!\n";
            }
            $feat{$feature} = $genes{$feature}{'start'};
        }

        my @sorted_features = sort {$feat{$a} <=> $feat{$b}}(keys %feat);

        foreach my $feature (@sorted_features){

            ## Gene features
            my $strand = $genes{$feature}{'strand'};
            my $gs = $genes{$feature}{'start'};
            my $ge = $genes{$feature}{'end'};
            my $lc = $genes{$feature}{'locus'};
            my $ftype = $genes{$feature}{'type'};

            if ($strand eq '+'){
                print GBFF "     gene            ".$gs.'..'.$ge."\n";
            }
            elsif ($strand eq '-'){
                print GBFF "     gene            ".'complement('.$gs.'..'.$ge.")\n";
            }
            print GBFF "                     /locus_tag=\"$lc\"\n";

            ## RNA features
            foreach my $rchild (@{$genes{$feature}{'children'}}){

                my $rproduct = $rnas{$rchild}{'product'};
                if ($id_flag){
                    $rproduct = $feature;
                }

                my @cd = @{$rnas{$rchild}{'coordinates'}};

                my %hash;
                for my $cd (@cd){
                    my ($key,$value) = $cd =~ /^(\d+)\.\.(\d+)$/;
                    $hash{$key} = $value;
                }

                my @scd;
                for my $key (sort {$a <=> $b}(keys %hash)){
                    my $value = $key.'..'.$hash{$key};
                    push (@scd, $value);
                }

                if (scalar(@scd) == 1){
                    if ($strand eq '+'){
                        print GBFF "     $ftype            ".$scd[0]."\n";
                    }
                    else {
                        print GBFF "     $ftype            ".'complement('.$scd[0].")\n";
                    }
                }
                else{

                    if ($strand eq '+'){
                        print GBFF "     $ftype            ".'join(';
                    }
                    else {
                        print GBFF "     $ftype            ".'complement(join(';
                    }

                    for (0..$#scd - 1){
                        print GBFF "$scd[$_],";
                    }
                    print GBFF "$scd[-1]";

                    if ($strand eq '+'){
                        print GBFF ")\n";
                    }
                    else {
                        print GBFF "))\n";
                    }

                }
                
                print GBFF "                     /locus_tag=\"$lc\"\n";
                print GBFF "                     /product=\"$rproduct\"\n";

                ## Adding CDS info
                if ($ftype eq 'mRNA'){

                    my @cd_cds = @{$rnas{$rchild}{'coordinates_cds'}};

                    my %hash;
                    for my $cd (@cd_cds){
                        my ($key,$value) = $cd =~ /^(\d+)\.\.(\d+)$/;
                        $hash{$key} = $value;
                    }

                    my @scd;
                    for my $key (sort {$a <=> $b}(keys %hash)){
                        my $value = $key.'..'.$hash{$key};
                        push (@scd, $value);
                    }

                    if (scalar(@scd) == 1){
                        if ($strand eq '+'){
                            print GBFF "     CDS             ".$scd[0]."\n";
                        }
                        else {
                            print GBFF "     CDS             ".'complement('.$scd[0].")\n";
                        }
                    }
                    else{
                        if ($strand eq '+'){
                            print GBFF "     CDS             ".'join(';
                        }
                        else {
                            print GBFF "     CDS             ".'complement(join(';
                        }

                        for (0..$#scd - 1){
                            print GBFF "$scd[$_],";
                        }
                        print GBFF "$scd[-1]";
                        
                        if ($strand eq '+'){
                            print GBFF ")\n";
                        }
                        else {
                            print GBFF "))\n";
                        }

                    }
                    
                    print GBFF "                     /locus_tag=\"$lc\"\n";
                    print GBFF "                     /product=\"$rproduct\"\n";

                    my $mRNA_substring;
                    foreach my $coord (@scd){
                        my ($start,$end) = $coord =~ /^(\d+)\.\.(\d+)$/;
                        my $len = $end - $start + 1;
                        $mRNA_substring .= substr($sequences{$sequence}, $start - 1, $len);
                    }

                    if ($strand eq '-'){
                        reverse_complement($mRNA_substring);
                    }

                    my $CDS = translate($mRNA_substring);
                    my $CDS_len = length($CDS);

                    if ($CDS_len <= 43){
                        print GBFF "                     /translation=\"".$CDS."\"\n";
                    }
                    else{

                        my $lCDS = substr($CDS, 0, 44);
                        my $rCDS = substr($CDS, 44);

                        print GBFF "                     /translation=\"".$lCDS."\n";

                        my @SEQUENCE = unpack ("(A58)*", $rCDS);
                        while (my $seq = shift@SEQUENCE){
                            if (length($seq) == 58){
                                print GBFF "                     $seq\n";
                            }
                            else {
                                print GBFF "                     $seq\"\n";
                            }
                        }

                    }

                }

            }

        }

        ## GBFF sequence
        print GBFF "ORIGIN      "."\n";
        my @sequence = unpack ("(A60)*", lc($sequences{$sequence}));
        my $seq_counter = 1;

        while (my $wide_60 = shift@sequence){

            my @wide_10 = unpack ("(A10)*", $wide_60);
            my $seq_counter_len = length($seq_counter);
            my $spacer_len = 9 - $seq_counter_len;
            my $spacer = ' ' x $spacer_len;

            print GBFF $spacer.$seq_counter;
            while (my $wide_10 = shift@wide_10){
                print GBFF ' '.$wide_10;
            }
            print GBFF "\n";

            $seq_counter += 60;

        }

        print GBFF '//'."\n";

    }

    ## Closing filehandles
    if ($gzip eq ':gzip'){
        binmode FASTA, ":gzip(none)";
    }

    close FASTA;
    close GBFF;

    ## Compress GFF files
    if ($gzip_flag){
        if ($verbose){
            print "Compressing ($gzip_tool): ".$outfile."\n";
        }
        system ("$gzip_tool $outfile");
    }

}


#########################################################################
### Subroutines
#########################################################################

sub reverse_complement {
    $_[0] = reverse($_[0]);
    $_[0] =~ tr/ATGCRYSWKMBDHVatgcryswkmbdhv/TACGYRWSMKVHDBtacgyrwsmkvhdb/;
}

sub translate {

    my $seq = uc($_[0]);
    my $protein;

    for (my $i = 0; $i < (length($seq) - 5); $i += 3){

        my $codon = substr($seq, $i, 3);

        if (exists $gcodes{$gc}{$codon}){
            $protein .=  $gcodes{$gc}{$codon};
        }
        else {
            $protein .= 'X';
        }
    }

    return $protein;

}

sub gcodes { ## NCBI Genetic codes
    %gcodes = (
        1 => { ## The Standard Code (transl_table=1)
            'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
            'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
            'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => '*',   
            'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',
            'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
            'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
            'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
            'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
            'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
            'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
            'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
            'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
            'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
            'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
            'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
            'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
        },
        2 => { ## The Vertebrate Mitochondrial Code (transl_table=2)
            'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',  
            'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
            'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => 'W',
            'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',
            'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
            'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
            'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
            'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
            'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
            'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
            'ATA' => 'M', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => '*',
            'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => '*',
            'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
            'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
            'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
            'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
        },
        3 => { ## The Yeast Mitochondrial Code (transl_table=3)
            'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
            'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
            'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => 'W',
            'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',
            'CTT' => 'T', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
            'CTC' => 'T', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
            'CTA' => 'T', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
            'CTG' => 'T', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
            'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
            'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
            'ATA' => 'M', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
            'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
            'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
            'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
            'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
            'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
        },
        4 => { ## The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code (transl_table=4)
            'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',  
            'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',  
            'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => 'W',  
            'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',  
            'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',  
            'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',  
            'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',  
            'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',  
            'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',  
            'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',  
            'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',  
            'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',  
            'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',  
            'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',  
            'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',  
            'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G' 
        },
        5 => { ## The Invertebrate Mitochondrial Code (transl_table=5)
            'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',  
            'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',  
            'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => 'W',  
            'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',  
            'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',  
            'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',  
            'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',  
            'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',  
            'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',  
            'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',  
            'ATA' => 'M', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'S',  
            'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'S',  
            'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',  
            'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',  
            'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',  
            'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'  
        },
        6 => { ## The Ciliate, Dasycladacean and Hexamita Nuclear Code (transl_table=6)
            'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C', 
            'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C', 
            'TTA' => 'L', 'TCA' => 'S', 'TAA' => 'Q', 'TGA' => '*', 
            'TTG' => 'L', 'TCG' => 'S', 'TAG' => 'Q', 'TGG' => 'W', 
            'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R', 
            'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R', 
            'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R', 
            'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R', 
            'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S', 
            'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S', 
            'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R', 
            'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R', 
            'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G', 
            'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G', 
            'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G', 
            'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G' 
        },
        9 => { ## The Echinoderm and Flatworm Mitochondrial Code (transl_table=9)
            'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',  
            'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',  
            'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => 'W',  
            'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',  
            'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',  
            'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',  
            'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',  
            'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',  
            'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',  
            'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',  
            'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'N', 'AGA' => 'S',  
            'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'S',  
            'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',  
            'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
            'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',  
            'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'  
        },
        10 => { ## The Euplotid Nuclear Code (transl_table=10)
            'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
            'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
            'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => 'C',
            'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',
            'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
            'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
            'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
            'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
            'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
            'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
            'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
            'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
            'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
            'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
            'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
            'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
        },
        11 => { ## The Bacterial, Archaeal and Plant Plastid Code (transl_table=11)
            'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',  
            'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',  
            'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => '*',  
            'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',  
            'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',  
            'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',  
            'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',  
            'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',  
            'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',  
            'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',  
            'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',  
            'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',  
            'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',  
            'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',  
            'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',  
            'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'  
        },
        12 => { ## The Alternative Yeast Nuclear Code (transl_table=12)
            'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
            'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
            'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => '*',
            'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',
            'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
            'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
            'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
            'CTG' => 'S', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
            'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
            'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
            'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
            'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
            'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
            'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
            'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
            'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
        },
        13 => { ## The Ascidian Mitochondrial Code (transl_table=13)
            'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
            'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
            'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => 'W',
            'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',
            'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
            'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
            'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
            'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
            'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
            'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
            'ATA' => 'M', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'G',
            'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'G',
            'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
            'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
            'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
            'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G',
        },
        14 => { ## The Alternative Flatworm Mitochondrial Code (transl_table=14)
            'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
            'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
            'TTA' => 'L', 'TCA' => 'S', 'TAA' => 'Y', 'TGA' => 'W',
            'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',
            'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
            'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
            'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
            'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
            'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
            'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
            'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'N', 'AGA' => 'S',
            'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'S',
            'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
            'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
            'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
            'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G',
        },
        16 => { ## Chlorophycean Mitochondrial Code (transl_table=16)
            'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
            'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
            'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => '*',
            'TTG' => 'L', 'TCG' => 'S', 'TAG' => 'L', 'TGG' => 'W',
            'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
            'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
            'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
            'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
            'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
            'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
            'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
            'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
            'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
            'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
            'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
            'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
        },
        21 => { ## Trematode Mitochondrial Code (transl_table=21)
            'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',  
            'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',  
            'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => 'W',  
            'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',  
            'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',  
            'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',  
            'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',  
            'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',  
            'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',  
            'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',  
            'ATA' => 'M', 'ACA' => 'T', 'AAA' => 'N', 'AGA' => 'S',  
            'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'S',  
            'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',  
            'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',  
            'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',  
            'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'  
        },
        22 => { ## Scenedesmus obliquus Mitochondrial Code (transl_table=22)
            'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',  
            'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',  
            'TTA' => 'L', 'TCA' => '*', 'TAA' => '*', 'TGA' => '*',  
            'TTG' => 'L', 'TCG' => 'S', 'TAG' => 'L', 'TGG' => 'W',  
            'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',  
            'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',  
            'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',  
            'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',  
            'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',  
            'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',  
            'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',  
            'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',  
            'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',  
            'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',  
            'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',  
            'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'  
        },
        23 => { ## Thraustochytrium Mitochondrial Code (transl_table=23)
            'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
            'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
            'TTA' => '*', 'TCA' => 'S', 'TAA' => '*', 'TGA' => '*',
            'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',
            'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
            'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
            'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
            'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
            'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
            'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
            'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
            'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
            'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
            'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
            'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
            'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
        },
        24 => { ## Pterobranchia Mitochondrial Code (transl_table=24)
            'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
            'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
            'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => 'W',
            'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',
            'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
            'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
            'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
            'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
            'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
            'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
            'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'S',
            'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'K',
            'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
            'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
            'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
            'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
        },
        25 => { ## Candidate Division SR1 and Gracilibacteria Code (transl_table=25)
            'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
            'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
            'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => 'G',
            'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',
            'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
            'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
            'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
            'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
            'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
            'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
            'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
            'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
            'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
            'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
            'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
            'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
        },
        26 => { ## Pachysolen tannophilus Nuclear Code (transl_table=26)
            'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
            'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
            'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => '*',
            'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',
            'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
            'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
            'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
            'CTG' => 'A', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
            'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
            'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
            'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
            'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
            'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
            'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
            'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
            'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
        },
        27 => { ## Karyorelict Nuclear (transl_table=27)
            'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
            'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
            'TTA' => 'L', 'TCA' => 'S', 'TAA' => 'Q', 'TGA' => 'W',
            'TTG' => 'L', 'TCG' => 'S', 'TAG' => 'Q', 'TGG' => 'W',
            'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
            'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
            'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
            'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
            'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
            'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
            'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
            'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
            'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
            'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
            'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
            'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
        },
        28 => { ## Condylostoma Nuclear (transl_table=28)
            'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
            'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
            'TTA' => 'L', 'TCA' => 'S', 'TAA' => 'Q', 'TGA' => 'W',
            'TTG' => 'L', 'TCG' => 'S', 'TAG' => 'Q', 'TGG' => 'W',
            'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
            'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
            'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
            'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
            'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
            'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
            'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
            'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
            'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
            'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
            'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
            'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
        },
        29 => { ## Mesodinium Nuclear (transl_table=29)
            'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
            'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
            'TTA' => 'L', 'TCA' => 'S', 'TAA' => 'Y', 'TGA' => '*',
            'TTG' => 'L', 'TCG' => 'S', 'TAG' => 'Y', 'TGG' => 'W',
            'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
            'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
            'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
            'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
            'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
            'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
            'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
            'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
            'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
            'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
            'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
            'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
        },
        30 => { ## Peritrich Nuclear (transl_table=30)
            'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
            'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
            'TTA' => 'L', 'TCA' => 'S', 'TAA' => 'E', 'TGA' => '*',
            'TTG' => 'L', 'TCG' => 'S', 'TAG' => 'E', 'TGG' => 'W',
            'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
            'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
            'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
            'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
            'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
            'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
            'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
            'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
            'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
            'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
            'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
            'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
        },
        31 => { ## Blastocrithidia Nuclear (transl_table=31)
            'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
            'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
            'TTA' => 'L', 'TCA' => 'S', 'TAA' => 'E', 'TGA' => 'W',
            'TTG' => 'L', 'TCG' => 'S', 'TAG' => 'E', 'TGG' => 'W',
            'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
            'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
            'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
            'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
            'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
            'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
            'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
            'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
            'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
            'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
            'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
            'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
        },
    );
}