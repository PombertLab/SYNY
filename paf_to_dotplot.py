#!/usr/bin/env python3
## Pombert lab, 2024
version = '0.1d'
updated = '2024-03-17'
name = 'paf_to_dotplot.py'

import sys
import os
import matplotlib.pyplot as plt
import argparse
import seaborn as sns

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
-o (--outdir)   Output directory [Default: ./]
-u (--unit)     Units divider [Default: 1e3]
-h (--height)   Figure height in inches [Default: 10.8]
-w (--width)    Figure width in inches [Default: 19.2]
-c (--color)    Color [Default: blue]
-n (--noticks)  Turn off ticks on x and y axes
-a (--palette)  Use a color palette (e.g. Spectral) instead
                of a monochrome plot
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
cmd.add_argument("-u", "--unit", default='1e3')
cmd.add_argument("-h", "--height", default=10.8)
cmd.add_argument("-w", "--width", default=19.2)
cmd.add_argument("-c", "--color", default='blue')
cmd.add_argument("-a", "--palette")
cmd.add_argument("-n", "--noticks", action='store_true')
args = cmd.parse_args()

paf_files = args.paf
outdir = args.outdir
unit = args.unit
height = args.height
width = args.width
ccolor = args.color
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

            # Alignments includes gaps: no. of bases
            # in query/subject often differs =>
            # creating a slope so that numbers of
            # x and y variables are the same for plotting
            q_span = q_end - q_start + 1
            s_span = s_end - s_start + 1
            slope = s_span / q_span

            y_start = None

            if orientation == '+':
                y_start = s_start
            else:
                y_start = s_end

            for x in range(q_start, q_end, 1):
                dataframe[query][subject][x] = y_start
                if orientation == '+':
                    y_start = y_start + slope
                else:
                    y_start = y_start - slope

    ##### Plotting
    x_axes_total = int(len(query_len_dict))
    y_axes_total = int(len(subject_len_dict))
    subplots_total = x_axes_total * y_axes_total

    # Setting default image to widescreen by default
    plt.rcParams["figure.figsize"] = (width,height)
    palette = sns.color_palette(color_palette, len(query_len_dict))

    fig, axes = plt.subplots(y_axes_total, x_axes_total, sharex='col', sharey='row')
    basename = os.path.basename(paf)
    fig.suptitle(basename)

    if noticks:
        plt.setp(axes, xticks=[], yticks=[])

    plt.rcParams.update({'font.size': 8})

    print("\nWorking on", basename, '=> ', end='')
    print('total subplots:', subplots_total)

    ynum = 0
    xnum = 0
    cnum = 0

    ## Convert scientific notation to number
    divider = int(float(unit))

    for query in sorted(query_len_dict):

        for subject in reversed(sorted(subject_len_dict)):

            if color_palette:
                ccolor = palette[cnum]

            xmax = query_len_dict[query] / divider
            ymax = subject_len_dict[subject] / divider

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
            if subject in dataframe[query].keys():
                x1 = dataframe[query][subject].keys()
                y1 = dataframe[query][subject].values()
                axes[ynum,xnum].scatter([x / divider for x in x1], [y / divider for y in y1], s=1, color=ccolor)

            ynum += 1

        cnum += 1
        xnum += 1
        ynum = 0


    ## Writing to output file
    output = basename
    filename = outdir + '/' + output.rsplit('.', 1)[0] + f".{unit}" + '.png'
    print(f"Creating {filename}...")
    plt.savefig(filename)
