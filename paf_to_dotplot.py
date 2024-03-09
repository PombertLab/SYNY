#!/usr/bin/env python3
## Pombert lab, 2024
version = '0.1'
updated = '2024-03-09'
name = 'paf_to_dotplot.py'

import sys
import os
import matplotlib.pyplot as plt
import argparse

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
-h (--height)   Figure height in inches [Default: 10.8]
-w (--width)    Figure width in inches [Default: 19.2]
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
args = cmd.parse_args()

paf_files = args.paf
outdir = args.outdir
height = args.height
width = args.width

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

    ### Plotting

    x_axes_total = int(len(query_len_dict))
    y_axes_total = int(len(subject_len_dict))
    subplots_total = x_axes_total * y_axes_total

    # Setting default image to widescreen by default
    plt.rcParams["figure.figsize"] = (width,height)

    fig, axes = plt.subplots(y_axes_total, x_axes_total, sharex='col', sharey='row')
    basename = os.path.basename(paf)
    fig.suptitle(basename)

    print("\nWorking on", basename, '=> ', end='')
    print('total subplots:', subplots_total)

    plt.rcParams.update({'font.size': 8})

    ynum = 0
    xnum = 0
    for query in sorted(query_len_dict):

        for subject in reversed(sorted(subject_len_dict)):

            xmax = query_len_dict[query] / 1000
            ymax = subject_len_dict[subject] / 1000

            axes[ynum,xnum].set_xlim(1,xmax)
            axes[ynum,xnum].set_ylim(1,ymax)

            if (ynum + 1) == len(subject_len_dict):
                axes[ynum,xnum].set_xlabel(query, rotation=45, ha='right')

            if (xnum == 0):
                axes[ynum,xnum].set_ylabel(subject, rotation=45, ha='right')

            if subject in dataframe[query].keys():
                x1 = dataframe[query][subject].keys()
                y1 = dataframe[query][subject].values()
                axes[ynum,xnum].scatter([x / 1000 for x in x1], [y / 1000 for y in y1], s=1)

            ynum += 1
        
        xnum += 1
        ynum = 0


    output = basename
    filename = outdir + '/' + output.rsplit('.', 1)[0] + '.png'
    print(f"Creating {filename}...")
    plt.savefig(filename)
