#!/usr/bin/env python3
## Pombert lab, 2024

version = '0.5c'
updated = '2024-05-27'
name = 'paf_to_dotplot.py'

import sys
import os
import re
import matplotlib
import matplotlib.pyplot as plt
import argparse
import seaborn as sns
from multiprocessing import Pool, Value
matplotlib.use('agg')

################################################################################
## README
################################################################################

usage = f"""
NAME        {name}
VERSION     {version}
UPDATED     {updated}
SYNOPSIS    Create genome alignment dotplots from PAF files with matplotlib

REQS        matplotlib

COMMAND    {name} \\
            -p *.paf \\
            -o DOTPLOTS 

OPTIONS:
-p (--paf)      PAF file(s) to plot
-f (--fasta)    FASTA files used to generate PAF alignments
-o (--outdir)   Output directory [Default: ./]
-u (--unit)     Units divider [Default: 1e5]
-h (--height)   Figure height in inches [Default: 10.8]
-w (--width)    Figure width in inches [Default: 19.2]
-c (--color)    Color [Default: blue]
-n (--noticks)  Turn off ticks on x and y axes
-x (--affix)    Add affix to filename
-a (--palette)  Use a color palette (e.g. Spectral) instead
                of a monochrome plot
--wdis          Horizontal distance (width) between subplots [Default: 0.05]
--hdis          Vertical distance (height) between subplots [Default: 0.1]
--fontsize      Font size [Default: 8]
--threads       Number of threads to use [Default: 16]
--version       Show script version
"""

# Print custom message if argv is empty
if (len(sys.argv) <= 1):
    print(usage)
    exit(0)

################################################################################
## Create command lines switches
################################################################################

cmd = argparse.ArgumentParser(add_help=False)
cmd.add_argument("-p", "--paf", nargs='*')
cmd.add_argument("-f", "--fasta", nargs='*')
cmd.add_argument("-o", "--outdir", default='./')
cmd.add_argument("-u", "--unit", default='1e5')
cmd.add_argument("-h", "--height", default=10.8)
cmd.add_argument("-w", "--width", default=19.2)
cmd.add_argument("-c", "--color", default='blue')
cmd.add_argument("-x", "--affix")
cmd.add_argument("-a", "--palette")
cmd.add_argument("-n", "--noticks", action='store_true')
cmd.add_argument("--wdis", default=0.05)
cmd.add_argument("--hdis", default=0.1)
cmd.add_argument("--fontsize", default=8)
cmd.add_argument("--threads", default=16)
cmd.add_argument("--version", action='store_true')
args = cmd.parse_args()

paf_files = args.paf
fasta_files = args.fasta
outdir = args.outdir
unit = args.unit
height = args.height
width = args.width
ccolor = args.color
affix = args.affix
color_palette = args.palette
noticks = args.noticks
wdis = args.wdis
hdis = args.hdis
fontsize = int(args.fontsize)
threads = int(args.threads)
scversion = args.version

#########################################################################
### Version
#########################################################################

if scversion:
    print ("\b")
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

if fasta_files is not None:

    for fasta in fasta_files:

        basename = os.path.basename(fasta)
        x = re.search(r'(\S+)\.fasta', basename)
        if x:
            basename = x.group(1)
            len_dict[basename] = {}

        seqname = None
        with open(fasta) as f:

            for line in f:

                m = re.search(r'>(\S+)', line)

                if m:
                    seqname = m.group(1)
                    len_dict[basename][seqname] = 0
                else:
                    length = len(line)
                    len_dict[basename][seqname] += length

################################################################################
## Parsing and plotting PAF file(s) 
################################################################################

lsize = len(paf_files) * 2
counter = Value('i', 0)

def dotplot(paf):

    global counter
    basename = os.path.basename(paf)
    qfile = None
    sfile = None

    m = re.search(r'^(\w+)_vs_(\w+)', basename)
    if m:
        sfile = m.group(1)
        qfile = m.group(2)

    dataframe = {}
    query_len_dict = {}
    subject_len_dict = {}

    if fasta_files is not None:
        query_len_dict = len_dict[qfile]
        subject_len_dict = len_dict[sfile]

    with open(paf) as file:

        matchnum = 0

        for line in file:

            matchnum += 1

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

            if fasta_files is None:
                query_len_dict[query] = query_len
                subject_len_dict[subject] = subject_len

            if query not in dataframe:
                dataframe[query] = {}
            if subject not in dataframe[query]:
                dataframe[query][subject] = {}
            
            dataframe[query][subject][matchnum] = {}
            if orientation == '+':
                dataframe[query][subject][matchnum] = [(q_start,q_end),(s_start,s_end)]
            else:
                dataframe[query][subject][matchnum] = [(q_start,q_end),(s_end,s_start)]

    ##### Output file(s)
    output = basename

    global ccolor
    acolor = ccolor

    if color_palette:
        acolor = color_palette

    global affix
    if noticks:
        affix = 'noticks'

    suffix = ''
    if affix is not None:
        suffix = '.' + affix

    pngfile = pngdir + '/' + output.rsplit('.', 1)[0] + suffix + f".{unit}" + f".{width}x{height}" + f".{acolor}" + '.png'
    svgfile = svgdir + '/' + output.rsplit('.', 1)[0] + suffix + f".{unit}" + f".{width}x{height}" + f".{acolor}" + '.svg'

    ##### Plotting
    x_axes_total = int(len(query_len_dict))
    y_axes_total = int(len(subject_len_dict))
    subplots_total = x_axes_total * y_axes_total

    # Setting default image to widescreen by default
    plt.rcParams['svg.fonttype'] = 'none'
    plt.rcParams["figure.figsize"] = (width,height)
    plt.rcParams.update({'font.size': fontsize})

    # color palette
    palette = sns.color_palette(color_palette, len(query_len_dict))

    if noticks:
        plt.setp(axes, xticks=[], yticks=[])

    ## Convert scientific notation to number
    divider = int(float(unit))

    ## No subplot
    if subplots_total == 1:

        xlabel = list(query_len_dict.keys())[0]
        ylabel = list(subject_len_dict.keys())[0]
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        xmax = query_len_dict[xlabel] / divider
        ymax = subject_len_dict[ylabel] / divider

        plt.xlim(1,xmax)
        plt.ylim(1,ymax)

        if color_palette:
            ccolor = palette[0]

        for x in dataframe[xlabel][ylabel]:
            xmono = dataframe[xlabel][ylabel][x][0]
            ymono = dataframe[xlabel][ylabel][x][1]
            plt.plot([x / divider for x in xmono], [y / divider for y in ymono], color=ccolor)

    ## With subplots
    elif subplots_total > 1:

        fig, axes = plt.subplots(y_axes_total, x_axes_total, sharex='col', sharey='row')
        fig.suptitle(basename)

        ynum = 0
        xnum = 0
        cnum = 0
        znum = 0

        for query in sorted(query_len_dict):

            for subject in reversed(sorted(subject_len_dict)):

                if color_palette:
                    ccolor = palette[cnum]

                xmax = query_len_dict[query] / divider
                ymax = subject_len_dict[subject] / divider

                ## 1-dimensional array
                if ((x_axes_total == 1) or (y_axes_total == 1)):

                    axes[znum].set_xlim(1,xmax)
                    axes[znum].set_ylim(1,ymax)

                    if y_axes_total > 1:

                        axes[znum].set_ylabel(subject, rotation=45, ha='right')

                        ## x-axes labels and ticks
                        if (ynum + 1) == len(subject_len_dict):
                            axes[znum].set_xlabel(query, rotation=45, ha='right')

                        if (ynum + 1) < len(subject_len_dict):
                            axes[znum].get_xaxis().set_visible(False)

                    else:
                        axes[znum].set_xlabel(query, rotation=45, ha='right')
                        axes[znum].set_ylabel(subject, rotation=45, ha='right')

                        if xnum > 0:
                            axes[znum].get_yaxis().set_visible(False)

                ## multidimensional array
                else:

                    axes[ynum,xnum].set_xlim(1,xmax)
                    axes[ynum,xnum].set_ylim(1,ymax)

                    ## x-axes labels and ticks
                    if (ynum + 1) == len(subject_len_dict):
                        axes[ynum,xnum].set_xlabel(query, rotation=45, ha='right')

                    if (ynum + 1) < len(subject_len_dict):
                        axes[ynum,xnum].get_xaxis().set_visible(False)

                    # y-axes labels and ticks
                    if (xnum == 0):
                        axes[ynum,xnum].set_ylabel(subject, rotation=45, ha='right')

                    if (xnum > 0):
                        axes[ynum,xnum].get_yaxis().set_visible(False)

                # subplots data
                if query in dataframe:

                    if subject in dataframe[query]:

                        for x in dataframe[query][subject]:

                            x1 = dataframe[query][subject][x][0]
                            y1 = dataframe[query][subject][x][1]

                            # 1-dimensional array
                            if ((x_axes_total == 1) or (y_axes_total == 1)):
                                axes[znum].plot([x / divider for x in x1], [y / divider for y in y1], color=ccolor)

                            # multidimensional array
                            else:
                                axes[ynum,xnum].plot([x / divider for x in x1], [y / divider for y in y1], color=ccolor)

                        if ((x_axes_total == 1) or (y_axes_total == 1)):
                            znum += 1

                ynum += 1

            cnum += 1
            xnum += 1
            ynum = 0

        ## Reducing space between subplots
        plt.subplots_adjust(
            wspace=float(wdis),
            hspace=float(hdis)
        )

    ##### Writing to output file
    with counter.get_lock():
        counter.value += 1
    print(f"{counter.value} / {lsize} - plotting {pngfile}...")
    plt.savefig(pngfile)

    with counter.get_lock():
        counter.value += 1
    print(f"{counter.value} / {lsize} - plotting {svgfile}...")
    plt.savefig(svgfile)

    ## Close fig
    plt.clf()
    plt.cla()
    plt.close('all')

## Run
pool = Pool(threads)
pool.map(dotplot, paf_files)
