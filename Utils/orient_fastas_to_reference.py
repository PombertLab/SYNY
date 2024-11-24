#!/usr/bin/env python3
## Pombert lab, 2023

version = '0.3.5'
updated = '2024-11-24'
name = 'orient_fastas_to_reference.py'

from sys import argv
from argparse import ArgumentParser
from os.path import isdir,basename
from os import makedirs,system
from textwrap import wrap
from datetime import datetime

################################################################################
## README
################################################################################

usage = f"""
NAME        {name}
VERSION     {version}
UPDATED     {updated}
SYNOPSIS    Performs BLASTN searches against a reference assembly (in FASTA format),
            then reorders the sequences in the FASTA queries based on the matches
            found against the reference.

REQS        NCBI BLAST+ - https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

USAGE       {name} \\
              -f 50507/50507.processed.fasta \\
              -r REF/assembly.fasta

OPTIONS
-f (--fasta)         -f (--fasta)         FASTA files to reorder/reorient
-r (--ref)           Reference genome assembly
-o (--outdir)        Output directory [Default:'oriented_fastas']

-i (--min_pident)    Minimum percent identity to assign segment to reference [Default: 85%]
-a (--min_palign)    Minimum percent of the contig participating in alignment to assign segment to reference [Default: 20%]
-x (--max_overlp)    Maximum percent of alignment allowed to overlap a previous alignment to assign segment to reference [Default: 5%]
--version            Show script version
"""

# Print custom message if argv is empty
if (len(argv) <= 1):
    print(usage)
    exit()

################################################################################
## Create command lines switches
################################################################################

GetOptions = ArgumentParser()

GetOptions.add_argument("-f","--fasta",nargs='+',required=True)
GetOptions.add_argument("-r","--ref",required=True)
GetOptions.add_argument("-o","--outdir",default='oriented_fastas')
GetOptions.add_argument("-i","--min_pident",default=85,type=float,choices=[x for x in range(0,101)])
GetOptions.add_argument("-a","--min_palign",default=20,type=float,choices=[x for x in range(0,101)])
GetOptions.add_argument("-x","--max_overlp",default=5,type=float,choices=[x for x in range(0,101)])
GetOptions.add_argument("--version", action='store_true')

args = GetOptions.parse_args()

fastas = args.fasta
ref = args.ref
outdir = args.outdir
min_pident = args.min_pident
min_palign = args.min_palign
max_overlp = args.max_overlp
scversion = args.version

#########################################################################
### Version
#########################################################################

if scversion:
    print ("")
    print (f"Script:     {name}")
    print (f"Version:    {version}")
    print (f"Updated:    {updated}\n")
    exit(0)

#########################################################################
### Custom functions
#########################################################################

def reverse_complement(sequence, /):

    '''
    This function returns the reverse complement of a DNA sequence.
    Expects a DNA string as input (positional argument).
    '''

    bases = {
        'A':'T','a':'t',
        'T':'A','t':'a',
        'C':'G','c':'g',
        'G':'C','g':'c',
        'N':'N','n':'n'
    }

    seq = ""

    for base in sequence[::-1]:
        seq += bases[base]

    return seq

def BLASTN(*, query, subject, outfile):

    '''
    This function runs BLASTN searches with a modified output format 6:
    'qseqid sseqid length pident qstart qend qlen sstart send slen sstrand bitscore'
    Expecs query, subject and outfile as keywords.
    '''

    print(f"BLASTN: {query} vs. {subject}")

    system(f"""
        blastn \\
          -query {query} \\
          -subject {subject} \\
          -outfmt "6 qseqid sseqid length pident qstart qend qlen sstart send slen sstrand bitscore" \\
          -out {outfile} \\
          2>/dev/null
    """)
#   -max_target_seqs 1 \\

#########################################################################
### Get sequences from REFERENCE assembly
#########################################################################

reference_seqs = {}
locus = False

REF = open(ref,'r')
for line in REF:
    line = line.strip()
    if line[0] == ">":
        locus = line.split()[0][1:]
        reference_seqs[locus] = ''
    elif locus:
        reference_seqs[locus] += line
REF.close()

sseqids = [x for x in sorted(reference_seqs.keys())]

#########################################################################
### Iterate over provided FASTA files
#########################################################################

if not isdir(outdir):
    makedirs(outdir,mode=0o755)

for file in fastas:

    filename = basename(file).split(".")[0]

    print(f"\nWorking on {filename}")

    #####################################################################
    ### Retrieve sequences from FASTA
    #####################################################################

    sequences = {}
    locus = False

    FASTA = open(file,'r')
    for line in FASTA:

        line = line.strip()

        if line[0] == '>':
            locus = line[1:].split(" ")[0]
            sequences[locus] = ""

        elif locus:
            sequences[locus] += line

    FASTA.close()

    qseqids = [x for x in sorted(sequences.keys())]

    #####################################################################
    ### Create output directory
    #####################################################################

    temp_dir = f"{outdir}/{filename}"

    if not isdir(temp_dir):
        makedirs(temp_dir,mode=0o755)

    #####################################################################
    ### Perform BLAST of FASTA vs REFERENCE
    #####################################################################

    BLASTN(query=file,subject=ref,outfile=f"{temp_dir}/results.blastn.6")

    #####################################################################
    ### Retrieve BLAST hits between FASTA and REFERENCE
    #####################################################################

    hits = {}

    # Parsing BLAST output file
    BLAST = open(f"{temp_dir}/results.blastn.6",'r')

    for line in BLAST:

        data = line.strip().split("\t")

        qseqid = data[0]
        sseqid = data[1]
        length = int(data[2])
        pident = float(data[3])
        qstart = int(data[4])
        qend = int(data[5])
        qlen = int(data[6])
        sstart = int(data[7])
        send = int(data[8])
        slen = int(data[9])
        strand = data[10]
        bitscore = float(data[11])

        if qseqid not in hits.keys():
            hits[qseqid] = []

        hits[qseqid].append(
            {
                'bitscore':bitscore,
                'pident':pident,
                'sseqid':sseqid,
                'length':length,
                'qlen':qlen,
                'qstart':qstart,
                'qend':qend,
                'slen':slen,
                'sstart':sstart,
                'send':send,
                'strand':strand,
            }
        )

    BLAST.close()

    # Reorienting sequences based on BLAST matches
    orientations = {}
    ref_assignment = {x:-1 for x in qseqids}
    alignment = {}
    assigned_locations = {}

    for qseqid in sorted(hits.keys()):

        assigned_bps = {x:False for x in range(len(sequences[qseqid]))}

        # Sorting matches by bitscore, from best to worse
        for result in sorted(hits[qseqid],key=lambda x: x['bitscore'],reverse=True):

            sseqid = result['sseqid']
            qstart = result['qstart']
            qend = result['qend']
            qlen = result['qlen']
            sstart = result['sstart']
            send = result['send']
            slen = result['slen']
            pident = result['pident']
            length = result['length']
            strand = result['strand']

            occupied_count = 0

            for bp in range(qstart,qend):
                if assigned_bps[bp]:
                    occupied_count += 1

            if pident >= min_pident:

                if (occupied_count < (min_palign/100)*length) and ((abs(qstart-qend)+1) > (max_overlp/100)*qlen):

                    for x in range(qstart,qend):
                        assigned_bps[x] = True

                    if sseqid not in assigned_locations.keys():
                        assigned_locations[sseqid] = []

                    if qseqid not in orientations.keys():
                        orientations[qseqid] = strand
                        ref_assignment[qseqid] = sseqids.index(sseqid)

                    if strand != 'plus':
                        sstart,send = send,sstart

                    if strand != orientations[qseqid]:
                        print(f"Inversion between {qseqid} and {sseqid}")
                        qstart,qend = qend,qstart

                    if qseqid not in alignment.keys():
                        alignment[qseqid] = [qstart,qend]

                    assigned_locations[sseqid].append(
                        {
                            'qseqid':qseqid,
                            'qstart':qstart,
                            'qend':qend,
                            'sstart':sstart,
                            'send':send,
                            'pident':pident,
                            'qlen':qlen,
                            'slen':slen,
                            'length':length,
                        }
                    )

    # Writing sequences to output files
    ORIENT = open(f"{temp_dir}/{filename}.oriented.fasta",'w')
    UNMATCHED = open(f"{temp_dir}/{filename}.unmatched.fasta",'w')

    for qseqid in qseqids:

        seq = sequences[qseqid]

        if qseqid in orientations.keys():
            if orientations[qseqid] == 'minus':
                seq = reverse_complement(seq)

        seq = "\n".join(wrap(seq,60))

        if qseqid in orientations.keys():
            ORIENT.write(f">{qseqid}\n")
            ORIENT.write(f"{seq}\n")

        else:
            UNMATCHED.write(f">{qseqid}\n")
            UNMATCHED.write(f"{seq}\n")

    ORIENT.close()
    UNMATCHED.close()

    # Creating summary + Circos links file
    LINKS = open(f"{temp_dir}/links.txt",'w')
    REF_MAP = open(f"{temp_dir}/all.map",'w')

    REF_MAP.write("## >REFERENCE_HIT\tREFERENCE_LENGTH\n")
    REF_MAP.write("##  >>FASTA_HIT\tALIGN_TYPE\tPIDENT\tFASTA_HIT_START\tFASTA_HIT_END\tFRACTION_FASTA_ALIGNED\tPERCENTAGE_FASTA_ALIGNED")
    REF_MAP.write("\tREFERENCE_HIT_START\tREFERENCE_HIT_END\tFRACTION_REFERENCE_ALIGNED\tPERCENTAGE_REFERENCE_ALIGNED\n\n")

    for ref_hit in sorted(assigned_locations.keys()):

        REF_MAP.write(f">>{ref_hit}\t{len(reference_seqs[ref_hit])}\n")

        for assignment in sorted(assigned_locations[ref_hit],key = lambda x: x['sstart']):

            pident = assignment['pident']
            qseqid = assignment['qseqid']
            qstart = assignment['qstart']
            qend = assignment['qend']
            sstart = assignment['sstart']
            send = assignment['send']
            qlen = assignment['qlen']
            slen = assignment['slen']
            length = assignment['length']

            aligned_bases = abs(qend-qstart)+1
            aligned_percent = aligned_bases/qlen*100
            reference_covered = length/slen*100

            REF_MAP.write(f" >{qseqid}")

            if qstart == alignment[qseqid][0] and qend == alignment[qseqid][-1]:
                REF_MAP.write(f"\tPrimary")
                LINKS.write(f"chr{(sseqids.index(ref_hit))+1} {sstart} {send} con{(qseqids.index(qseqid))+1} {qstart} {qend} color=0,255,255,.25,z=0\n")

            else:
                REF_MAP.write(f"\tSecondary")
                LINKS.write(f"chr{(sseqids.index(ref_hit))+1} {sstart} {send} con{(qseqids.index(qseqid))+1} {qstart} {qend} color=255,0,255,.25,z=10\n")

            REF_MAP.write(f"\t{pident}%")
            REF_MAP.write(f"\t{qstart}")
            REF_MAP.write(f"\t{qend}")
            REF_MAP.write(f"\t{aligned_bases}/{qlen}")
            REF_MAP.write(f"\t{aligned_percent:.2f}%")
            REF_MAP.write(f"\t{sstart}")
            REF_MAP.write(f"\t{send}")
            REF_MAP.write(f"\t{length}/{slen}")
            REF_MAP.write(f"\t{reference_covered:.2f}%\n")

        REF_MAP.write("\n")

    REF_MAP.close()
    LINKS.close()

    # Creating Circos karyotype file
    KARYO = open(f"{temp_dir}/karyotype.txt",'w')

    KARYO.write(f"# reference karyotype\n")
    for index,key in enumerate(sorted(reference_seqs.keys())):
        KARYO.write(f"chr - chr{index+1} {key} 0 {len(reference_seqs[key])} chr1\n")

    KARYO.write(f"\n# assembly karyotype\n")
    for key in sorted(ref_assignment.keys(),key = lambda x: ref_assignment[x],reverse=True):
        index = qseqids.index(key)
        KARYO.write(f"chr - con{index+1} {key} 0 {len(sequences[key])} chr5\n")

### Completion
time = datetime.now()
print("\n", 'Process completed on: ', time.strftime("%d/%m/%Y %H:%M:%S"), "\n", sep='')