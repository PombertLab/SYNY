#!/usr/bin/env python3
## Pombert lab, 2024
version = '0.1a'
updated = '2024-03-16'
name = 'paf_to_barplot.py'

import sys
import os
import argparse
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
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
            # See https://www.practicalpythonfordatascience.com/ap_seaborn_palette
            # for a list of color palettes

COMMAND    {name} \\
            -p *.paf \\
            -o BARPLOTS 

OPTIONS:
-p (--paf)      PAF file(s) to plot
-o (--outdir)   Output directory [Default: ./]
-h (--height)   Figure height in inches [Default: 10.8]
-w (--width)    Figure width in inches [Default: 19.2]
-c (--palette)  Seaborn color palette [Default: Spectral]
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
cmd.add_argument("-o", "--outdir", default='./')
cmd.add_argument("-h", "--height", default=10.8)
cmd.add_argument("-w", "--width", default=19.2)
cmd.add_argument("-c", "--palette", default='Spectral')
cmd.add_argument("-n", "--noticks", action='store_true')
args = cmd.parse_args()

paf_files = args.paf
outdir = args.outdir
height = args.height
width = args.width
color_palette = args.palette
noticks = args.noticks

################################################################################
## Working on output directory
################################################################################

if os.path.isdir(outdir) == False:
    try:
        os.makedirs(outdir)
    except:
        sys.exit(f"Can't create directory {outdir}...")

################################################################################
## Parsing and plotting PAF file(s) 
################################################################################

for paf in paf_files:

    query_len_dict = {}
    subject_len_dict = {}
    dataframe = {}

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

            query_len_dict[query] = query_len
            subject_len_dict[subject] = subject_len

            if query not in dataframe:
                dataframe[query] = {}
            if subject not in dataframe[query]:
                dataframe[query][subject] = {}

            y_width = s_end - s_start + 1
            dataframe[query][subject][s_start] = y_width


    ### Plotting 
    palette = sns.color_palette(color_palette, len(query_len_dict))

    # Setting default image to widescreen by default
    plt.rcParams["figure.figsize"] = (width,height)
    plt.rcParams.update({'font.size': 8})


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

    ## 
    x = 0.5
    cnum = 0
    legend = []
    dlegend = {}
    for subject in reversed(sorted(subject_len_dict.keys())):

        for query in dataframe.keys():

            ccolor=palette[cnum]
            xlegend = mpatches.Patch(color=ccolor, label=query)
            if query not in dlegend:
                dlegend[query] = {}
                legend.append(xlegend)

            cnum += 1

            if subject in dataframe[query].keys():

                for start in dataframe[query][subject].keys():
                    end = dataframe[query][subject][start]
                    ax.add_patch(Rectangle((start, x), end, 1, color=ccolor))

        x += 2
        cnum = 0

    ## Adding labels to y-axis
    ax.set_yticks(ydata)
    ax.set_yticks(ax.get_yticks()[:-1])
    ax.set_yticklabels(ylabels)


    ## Writing to output file
    basename = os.path.basename(paf)
    output = basename
    fig.suptitle(basename)

    plt.legend(handles=legend, loc='center left', bbox_to_anchor=(1, 0.5))
    filename = outdir + '/' + output.rsplit('.', 1)[0] + '.barplot.png'
    svg = outdir + '/' + output.rsplit('.', 1)[0] + '.barplot.svg'

    print(f"Creating {filename}...")
    plt.savefig(filename)

    print(f"Creating {svg}...")
    plt.savefig(svg)