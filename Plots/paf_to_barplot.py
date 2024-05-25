#!/usr/bin/env python3
## Pombert lab, 2024
version = '0.4'
updated = '2024-05-25'
name = 'paf_to_barplot.py'

import sys
import os
import re
import argparse
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
from multiprocessing import Pool, Value
import seaborn as sns

################################################################################
## README
################################################################################

usage = f"""
NAME        {name}
VERSION     {version}
UPDATED     {updated}
SYNOPSIS    Create genome alignment dotplots from PAF files with matplotlib

REQS        matplotlib, seaborn

COMMAND    {name} \\
            -p *.paf \\
            -o BARPLOTS 

OPTIONS:
-p (--paf)      PAF file(s) to plot
-f (--fasta)    FASTA files used to generate PAF alignments
-o (--outdir)   Output directory [Default: ./]
-h (--height)   Figure height in inches [Default: 10.8]
-w (--width)    Figure width in inches [Default: 19.2]
-a (--affix)    Add affix to filename
-c (--palette)  Seaborn color palette [Default: Spectral]
                # See https://www.practicalpythonfordatascience.com/ap_seaborn_palette
                # for a list of color palettes
-m (--mono)     Use a single specified mochochrome color for all chromosomes instead
                of a color palette
--clusters      Color matches by clusters instead; colors are not related between contigs
--catenate      Concatenate queries against the same reference in a single barplot
--fontsize      Font size [Default: 8]
--threads       Number of threads to use [Default: 16]
"""

# Print custom message if argv is empty
if (len(sys.argv) <= 1):
    print(usage)
    sys.exit()

################################################################################
## Create command lines switches
################################################################################

cmd = argparse.ArgumentParser(add_help=False)
cmd.add_argument("-p", "--paf", nargs='*')
cmd.add_argument("-f", "--fasta", nargs='*')
cmd.add_argument("-o", "--outdir", default='./')
cmd.add_argument("-h", "--height", default=10.8)
cmd.add_argument("-w", "--width", default=19.2)
cmd.add_argument("-a", "--affix")
cmd.add_argument("-c", "--palette", default='Spectral')
cmd.add_argument("-n", "--noticks", action='store_true')
cmd.add_argument("-m", "--mono")
cmd.add_argument("--clusters", action='store_true')
cmd.add_argument("--catenate", action='store_true')
cmd.add_argument("--fontsize", default=8)
cmd.add_argument("--threads", default=16)
args = cmd.parse_args()

paf_files = args.paf
fasta_files = args.fasta
outdir = args.outdir
height = args.height
width = args.width
affix = args.affix
color_palette = args.palette
noticks = args.noticks
monochrome = args.mono
clusters = args.clusters #
catenate = args.catenate
fontsize = int(args.fontsize)
threads = int(args.threads)

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
cat_dict = {}

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
                    cat_dict[seqname] = 0
                else:
                    length = len(line)
                    len_dict[basename][seqname] += length
                    cat_dict[seqname] += length

################################################################################
## Parsing and plotting PAF file(s) 
################################################################################

lsize = len(paf_files) * 2
counter = Value('i', 0)

def barplot(paf):

    basename = os.path.basename(paf)
    qfile = None
    sfile = None

    global counter

    m = re.search(r'^(\w+)_vs_(\w+)', basename)
    if m:
        sfile = m.group(1)
        qfile = m.group(2)

    dataframe = {}
    query_len_dict = {}
    subject_len_dict = {}

    if fasta_files is not None:
        query_len_dict = len_dict[qfile]
        if catenate:
            for key in len_dict:
                if key != qfile:
                    subject_len_dict.update(len_dict[key])
        else:
            subject_len_dict = len_dict[sfile]

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

            if fasta_files is None:
                query_len_dict[query] = query_len
                subject_len_dict[subject] = subject_len

            if query not in dataframe:
                dataframe[query] = {}
            if subject not in dataframe[query]:
                dataframe[query][subject] = {}

            y_width = s_end - s_start + 1
            dataframe[query][subject][s_start] = y_width

    ##### Plotting data
    palette = sns.color_palette(color_palette, len(query_len_dict))

    # Setting default image to widescreen by default
    plt.rcParams['svg.fonttype'] = 'none'
    plt.rcParams["figure.figsize"] = (width,height)
    plt.rcParams.update({'font.size': fontsize})

    # Defining Matplotlib figure and axis
    fig, ax = plt.subplots()

    if noticks:
        plt.setp(ax, xticks=[], yticks=[])

    # creating a blank histogram
    ax.hist([0, 0],[0, 0], orientation='horizontal', color='white')

    ## Creating bar plot edges
    x = 0.5
    ydata = [1]
    ylabels = []
    for key in reversed(sorted(subject_len_dict.keys())):
        ylabels.append(key)
        ax.add_patch(Rectangle((1, x), subject_len_dict[key], 1, edgecolor='black', fc='none'))
        x += 2
        ydata.append(x + 0.5)

    ## Filling bar plots
    x = 0.5
    cnum = 0
    legend = []
    dlegend = {}

    for subject in reversed(sorted(subject_len_dict.keys())):

        for query in sorted(query_len_dict.keys()):

            if not clusters:
                ccolor = palette[cnum]
                xlegend = mpatches.Patch(color=ccolor, label=query)
                if query not in dlegend:
                    dlegend[query] = {}
                    legend.append(xlegend)

            if monochrome:
                ccolor = monochrome

            cnum += 1

            if query in dataframe:
                xnum = 0
                if subject in dataframe[query]:
                    for start in dataframe[query][subject]:
                        end = dataframe[query][subject][start]
                        if clusters:
                            palette = sns.color_palette(color_palette, len(dataframe[query][subject]))
                            ccolor = palette[xnum]
                            ax.add_patch(Rectangle((start, x), end, 1, color=ccolor))
                        else:
                            ax.add_patch(Rectangle((start, x), end, 1, color=ccolor))
                        xnum += 1

        x += 2
        cnum = 0

    ## Adding labels to y-axis
    ax.set_yticks(ydata)
    ax.set_yticks(ax.get_yticks()[:-1])
    ax.set_yticklabels(ylabels)

    ## Hiding border frame
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # ax.spines['bottom'].set_visible(False)

    ## Writing to output file
    basename = os.path.basename(paf)
    output = basename
    fig.suptitle(basename)

    if monochrome or clusters:
        next
    else:
        plt.legend(handles=legend, loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)


    affix_color = color_palette
    if monochrome:
        affix_color = monochrome
    
    if clusters:
        affix_color = 'clus'

    suffix = ''
    if affix is not None:
        suffix = '.' + affix

    png = pngdir + '/' + output.rsplit('.', 1)[0] + suffix + '.barplot.' + f"{width}x{height}." + f"{affix_color}" + '.png'
    svg = svgdir + '/' + output.rsplit('.', 1)[0] + suffix + '.barplot.' + f"{width}x{height}." + f"{affix_color}" + '.svg'

    with counter.get_lock():
        counter.value += 1
    print(f"{counter.value} / {lsize} - plotting {png}...")
    if monochrome or clusters:
        plt.savefig(png)
    else:
        plt.savefig(png, bbox_inches='tight')

    with counter.get_lock():
        counter.value += 1
    print(f"{counter.value} / {lsize} - plotting {svg}...")
    if monochrome or clusters:
        plt.savefig(svg)
    else:
        plt.savefig(svg, bbox_inches='tight')

    ## Close fig
    plt.clf()
    plt.cla()
    plt.close('all')

## Run
pool = Pool(threads)
pool.map(barplot, paf_files)