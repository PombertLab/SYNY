#!/usr/bin/env python3
## Pombert lab, 2024

version = '0.3'
updated = '2024-03-25'
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
-u (--unit)     Units divider [Default: 1e5]
-h (--height)   Figure height in inches [Default: 10.8]
-w (--width)    Figure width in inches [Default: 19.2]
-c (--color)    Color [Default: blue]
-n (--noticks)  Turn off ticks on x and y axes
-a (--palette)  Use a color palette (e.g. Spectral) instead
                of a monochrome plot
--wdis          Horizontal distance (width) between subplots [Default: 0.05]
--hdis          Vertical distance (height) between subplots [Default: 0.1]
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
cmd.add_argument("-u", "--unit", default='1e5')
cmd.add_argument("-h", "--height", default=10.8)
cmd.add_argument("-w", "--width", default=19.2)
cmd.add_argument("-c", "--color", default='blue')
cmd.add_argument("-a", "--palette")
cmd.add_argument("-n", "--noticks", action='store_true')
cmd.add_argument("--wdis", default=0.05)
cmd.add_argument("--hdis", default=0.1)
args = cmd.parse_args()

paf_files = args.paf
outdir = args.outdir
unit = args.unit
height = args.height
width = args.width
ccolor = args.color
color_palette = args.palette
noticks = args.noticks
wdis = args.wdis
hdis = args.hdis

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

        basename = os.path.basename(paf)
        print(f"\nWorking on {basename}", end='')

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

            if orientation == '-':
                slope = -abs(slope)
                s_start = s_end

            dataframe[query][subject][q_start] = [q_end,s_start,slope]

    ##### Plotting
    x_axes_total = int(len(query_len_dict))
    y_axes_total = int(len(subject_len_dict))
    subplots_total = x_axes_total * y_axes_total
    print( '  => ', 'total subplots:', subplots_total)

    # Setting default image to widescreen by default
    plt.rcParams["figure.figsize"] = (width,height)
    plt.rcParams.update({'font.size': 8})

    # color palette
    palette = sns.color_palette(color_palette, len(query_len_dict))

    if noticks:
        plt.setp(axes, xticks=[], yticks=[])

    ## Convert scientific notation to number
    divider = int(float(unit))

    fig = None
    axes = None

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

        xmono = []
        ymono = []

        for q_start in dataframe[xlabel][ylabel].keys():

            q_end = dataframe[xlabel][ylabel][q_start][0]
            s_start = dataframe[xlabel][ylabel][q_start][1]
            slope = dataframe[xlabel][ylabel][q_start][2]

            y = s_start
            for x in range(q_start,q_end):
                y += slope
                xmono.append(x)
                ymono.append(y)

        fig = plt.scatter([x / divider for x in xmono], [y / divider for y in ymono], s=1, color=ccolor)

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
                if subject in dataframe[query].keys():

                    x1 = []
                    y1 = []

                    for q_start in dataframe[query][subject].keys():

                        q_end = dataframe[query][subject][q_start][0]
                        s_start = dataframe[query][subject][q_start][1]
                        slope = dataframe[query][subject][q_start][2]

                        y = s_start
                        for x in range(q_start,q_end):
                            y += slope
                            x1.append(x)
                            y1.append(y)

                    # 1-dimensional array
                    if ((x_axes_total == 1) or (y_axes_total == 1)):
                        axes[znum].scatter([x / divider for x in x1], [y / divider for y in y1], s=1, color=ccolor)
                        znum += 1

                    # multidimensional array
                    else:
                        axes[ynum,xnum].scatter([x / divider for x in x1], [y / divider for y in y1], s=1, color=ccolor)

                    ## To help reduce memory footprint
                    del dataframe[query][subject]

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
    output = basename

    acolor = ccolor
    affix = unit

    if color_palette:
        acolor = color_palette

    if noticks:
        affix = 'noticks'

    filename = outdir + '/' + output.rsplit('.', 1)[0] + f".{affix}" + f".{width}x{height}" + f".{acolor}" + '.png'
    print(f"Creating {filename}...")
    plt.savefig(filename)

    ## Close fig
    plt.clf()
    plt.close()