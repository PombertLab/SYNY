#!/usr/bin/env python3
## Pombert lab, 2024
version = '0.1b'
updated = '2024-05-27'
name = 'linear_maps.py'

import sys
import os
import re
import argparse
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from multiprocessing import Pool, Value
import seaborn as sns

################################################################################
## README
################################################################################

usage = f"""
NAME        {name}
VERSION     {version}
UPDATED     {updated}
SYNOPSIS    Create linear maps from PAF files with matplotlib

REQS        matplotlib

COMMAND    {name} \\
            -f *.fasta \\
            -p *.paf \\
            -o LINEMAPS

OPTIONS:
-f (--fasta)        FASTA files used to generate PAF alignments
-p (--paf)          PAF file(s) to parse
-o (--outdir)       Output directory [Default: ./]
-c (--rpalette)     Color palette to use for ref genome [Default: Spectral]
-x (--xpalette)     Color palette to use for target genome [Default: Blues]
-r (--rotation)     Contig name rotation [Default: 90]
-h (--height)       Set figure height [Default: 5]
-w (--width)        Set figure width [Default: 20]
--fontsize          Font size [Default: 8]
--threads           Number of threads to use [Default: 16]
--version           Show script version
"""

# Print custom message if argv is empty
if (len(sys.argv) <= 1):
    print(usage)
    exit(0)

################################################################################
## Create command lines switches
################################################################################

cmd = argparse.ArgumentParser(add_help=False)
cmd.add_argument("-f", "--fasta", nargs='*')
cmd.add_argument("-p", "--paf", nargs='*')
cmd.add_argument("-o", "--outdir", default='./')
cmd.add_argument("-c", "--rpalette", default='Spectral')
cmd.add_argument("-x", "--xpalette", default='Blues')
cmd.add_argument("-r", "--rotation", default=90)
cmd.add_argument("-h", "--height", default=5)
cmd.add_argument("-w", "--width", default=20)
cmd.add_argument("--fontsize", default=8)
cmd.add_argument("--threads", default=16)
cmd.add_argument("--version", action='store_true')
args = cmd.parse_args()

fasta_files = args.fasta
paf_files = args.paf
outdir = args.outdir
rpalette = args.rpalette
xpalette = args.xpalette
rotation = int(args.rotation)
height = int(args.height)
width = int(args.width)
fontsize = int(args.fontsize)
threads = int(args.threads)
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

################################################################################
## Working on output directory
################################################################################

svgdir = outdir + '/SVG'
pngdir = outdir + '/PNG'

subdirs = [outdir,svgdir,pngdir]

for dir in subdirs:
    if os.path.isdir(dir) == False:
        try:
            os.makedirs(dir)
        except:
            sys.exit(f"Can't create directory {dir}...")

################################################################################
## Parsing FASTA file(s) 
################################################################################

len_dict = {}
cumul_seq_len_dict = {}
tot_genome_len = {}

if fasta_files is not None:

    for fasta in fasta_files:

        genome_len = 0
        cumul_seq_length = 0

        basename = os.path.basename(fasta)
        x = re.search(r'(\S+)\.fasta', basename)
        if x:
            basename = x.group(1)
            len_dict[basename] = {}
            cumul_seq_len_dict[basename] = {}

        seqname = None
        with open(fasta) as f:

            for line in f:

                m = re.search(r'>(\S+)', line)

                if m:
                    seqname = m.group(1)
                    len_dict[basename][seqname] = 0
                    cumul_seq_length = genome_len
                    cumul_seq_len_dict[basename][seqname] = genome_len
                else:
                    length = len(line)
                    len_dict[basename][seqname] += length
                    genome_len += length

        tot_genome_len[basename] = genome_len


################################################################################
## Parsing PAF file(s) 
################################################################################

lsize = len(paf_files) * 2
counter = Value('i', 0)

def linemap(paf):

    dataframe = {}
    points = {}

    basename = os.path.basename(paf)
    qfile = None
    sfile = None
    ref = None
    target = None
    global counter

    m = re.search(r'^(\w+)_vs_(\w+)', basename)
    if m:
        sfile = m.group(1)
        qfile = m.group(2)
        target = m.group(1)
        ref = m.group(2)

    with open(paf) as file:

        for line in file:

            # PAF data structure:
            # https://github.com/lh3/miniasm/blob/master/PAF.md

            data = line.split("\t")
            query = data[0]
            query_len = int(data[1])
            q_start = int(data[2])
            q_end = int(data[3])
            orientation = data[4]
            subject = data[5]
            subject_len = int(data[6])
            s_start = int(data[7])
            s_end = int(data[8])

            if query not in dataframe:
                dataframe[query] = {}
            if subject not in dataframe[query]:
                dataframe[query][subject] = {}
            if query not in points:
                points[query] = []

            sstart = cumul_seq_len_dict[qfile][query] + q_start
            send = cumul_seq_len_dict[qfile][query] + q_end
            xstart = cumul_seq_len_dict[sfile][subject] + s_end
            xend  = cumul_seq_len_dict[sfile][subject] + s_start
            if orientation == '-':
                xstart = cumul_seq_len_dict[sfile][subject] + s_start
                xend = cumul_seq_len_dict[sfile][subject] + s_end

            ## To create polygon edges for ribbons between matches:
            ## If same orientation, x1 - x2, then y2 - y1
            ## If reversed, x1 - x2, then y1 - y2
            point = [
                [xstart, 0], ## x1 + 0 : x edge 1
                [xend, 0], ## x2 + 0 : x edge 2
                [sstart, 1], ## y2 + 1 : y edge 1 
                [send, 1], ## y1 + 1 : y edge 2
            ]
            points[query].append(point)


    ########## Plotting data

    # Setting default image settings
    plt.rcParams['svg.fonttype'] = 'none'
    plt.rcParams["figure.figsize"] = (width,height)
    plt.rcParams.update({'font.size': fontsize})

    ##### Karyogram
    pairwise = [ref, target]

    # Setting max length for x-axis
    maxlen = 0
    for x in pairwise:
        xlen = tot_genome_len[x]
        if xlen > maxlen:
            maxlen = xlen

    ## Defining karyotypes
    genomes = {}
    for genome in pairwise:

        genomes[genome] = {}

        cumulative_length = 0
        if genome == ref:
            palette = sns.color_palette(rpalette, len(len_dict[genome]))
        else:
            palette = sns.color_palette(xpalette, len(len_dict[genome]))
        cnum = 0

        for contig in len_dict[genome]:

            genomes[genome][contig] = {}

            cwidth = len_dict[genome][contig] - 1
            genomes[genome][contig]['start'] = cumulative_length
            genomes[genome][contig]['end'] = cumulative_length + cwidth
            genomes[genome][contig]['color'] = palette[cnum]
            cumulative_length += cwidth + 1
            cnum += 1

    ## Creating 3 subplots; 2 karyotypes plus linsk (polygons) sandwiched in-between 
    y_axes_total = 3
    x_axes_total = 1
    fig, ax = plt.subplots(y_axes_total, x_axes_total, sharex='row', figsize=(width,height))

    def karyogram(bands, xnum):

        for band in bands:
            start, end, pcolor = band
            bwidth = end - start + 1
            ax[xnum].set_ylim(0, 1)
            ax[xnum].add_patch(plt.Rectangle((start, 0), bwidth, 1, color=pcolor))

        ax[xnum].set_xlim(0, maxlen)

        ## Hiding border frames
        ax[xnum].spines['top'].set_visible(False)
        ax[xnum].spines['left'].set_visible(False)
        ax[xnum].spines['right'].set_visible(False)
        ax[xnum].spines['bottom'].set_visible(False)

    ###### Drawing karyograms
    xnum = 0
    for genome in pairwise:

        bands = []
        ticks = []
        xlabels = list(genomes[genome])

        for contigs in genomes[genome]:

            start = genomes[genome][contigs]['start']
            end = genomes[genome][contigs]['end']
            color = genomes[genome][contigs]['color']

            coordinates = []
            coordinates.append(start)
            coordinates.append(end)
            coordinates.append(color)
            bands.append(coordinates)

            xwidth = end - start + 1
            midpoint = start + int(xwidth/2)
            ticks.append(midpoint)

        if xnum == 0:
            ax[xnum].tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)

        ax[xnum].set_xticks(ticks, labels=xlabels, rotation=rotation)
        ax[xnum].set_yticks([])
        ax[xnum].set_ylabel(genome)
        karyogram(bands, xnum)
        xnum += 2


    ###### Drawing links (polygons) between karyograms
    ax[1].set_xlim(0, maxlen)

    pnum = 0
    for query in len_dict[ref]:

        palette = sns.color_palette(rpalette, len(list(len_dict[ref])))
        pcolor = palette[pnum]
        pnum += 1

        if query not in points:
            next
        else:
            for point in points[query]:
                polygon= plt.Polygon(point, fill=True, edgecolor=pcolor, facecolor=pcolor)
                ax[1].add_patch(polygon)

    ax[1].axis('off')

    plt.xlim(0, maxlen)
    plt.xticks(horizontalalignment='center')

    ## Creating Output files
    prefix = None
    m = re.search(r'^(.*?)\.paf', basename)
    if m:
        prefix = m.group(1)

    pngfile = pngdir + '/' + prefix + f".linemap" + f".{width}x{height}" + f".{rpalette}" + '.png'
    svgfile = svgdir + '/' + prefix + f".linemap" + f".{width}x{height}" + f".{rpalette}" + '.svg'

    with counter.get_lock():
        counter.value += 1
    print(f"{counter.value} / {lsize} - plotting {pngfile}...")
    plt.savefig(pngfile, bbox_inches='tight')

    with counter.get_lock():
        counter.value += 1
    print(f"{counter.value} / {lsize} - plotting {svgfile}...")
    plt.savefig(svgfile, bbox_inches='tight')

    ## Close fig
    plt.clf()
    plt.cla()
    plt.close('all')

## Run
pool = Pool(threads)
pool.map(linemap, paf_files)