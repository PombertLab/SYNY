#!/usr/bin/env python3
## Pombert lab, 2024

name = 'protein_cluster_hm.py'
version = '0.1a'
updated = '2024-03-27'

import sys
import os
import matplotlib
import matplotlib.pyplot as plt
import argparse
import re
import numpy as np
import seaborn as sns

################################################################################
## README
################################################################################

usage = f"""
NAME        {name}
VERSION     {version}
UPDATED     {updated}
SYNOPSIS    Create heatmaps from SYNY summary files with matplotlib

REQS        matplotlib

COMMAND    {name} \\
            -t clusters_summary_table.tsv \\
            -o HEATMAPS 

OPTIONS:
-t (--tsv)      clusters_summary_table.tsv
-o (--outdir)   Output directory [Default: ./HEATMAPS]
-p (--palette)  Color palette [Default: husl]
-c (--color)    Digit colors [Default: crest]
-h (--height)   Figure height in inches [Default: 10]
-w (--width)    Figure width in inches [Default: 10]
"""

# Print custom message if argv is empty
if (len(sys.argv) <= 1):
    print(usage)
    sys.exit()

################################################################################
## Create command lines switches
################################################################################

cmd = argparse.ArgumentParser(add_help=False)
cmd.add_argument("-t", "--tsv")
cmd.add_argument("-o", "--outdir", default='./HEATMAPS')
cmd.add_argument("-h", "--height", default=10)
cmd.add_argument("-w", "--width", default=10)
cmd.add_argument("-p", "--palette", default='crest')
cmd.add_argument("-c", "--color", default='white')
args = cmd.parse_args()

tsv_file = args.tsv
outdir = args.outdir
height = args.height
width = args.width
color_palette = args.palette
num_color = args.color

################################################################################
## Working on output directory
################################################################################

if os.path.isdir(outdir) == False:
    try:
        os.makedirs(outdir)
    except:
        sys.exit(f"Can't create directory {outdir}...")

################################################################################
## Parsing summary table 
################################################################################

dict = {}

with open (tsv_file) as f:

    for line in f:

        if '### Query' in line:
            next

        else:
            data = line.split("\t")
            comparison = data[0]
            gap = data[2]
            percent = data[4]

            if gap not in dict:
                dict[gap] = {}

            m = re.search(r'^(\w+)_vs_(\w+)$', comparison)
            if m:
                query = m.group(1)
                subject = m.group(2)

                if query not in dict[gap]:
                    dict[gap][query] = {}

                if subject not in dict[gap][query]:
                    dict[gap][query][subject] = {}

                ## Adding query against itself since skipped in comparisosn
                if query not in dict[gap][query]:
                    dict[gap][query][query] = 100

                dict[gap][query][subject] = percent

## Setting value as zero for missing matches in PAF
for gap in dict.keys():
    for query in dict[gap].keys():
        for subject in dict[gap].keys():
            if subject not in dict[gap][query]:
                dict[gap][query][subject] = 0


################################################################################
## Plotting summary table 
################################################################################

plt.rcParams["figure.figsize"] = (width,height)

for gap in dict.keys():

    ## Creating labels from keys
    labels = sorted(dict[gap].keys())

    ## Creating numpy array
    array = []
    for query in sorted(dict[gap].keys()):
        list = []
        for subject in sorted(dict[gap][query]):
            list.append(float(dict[gap][query][subject]))
        array.append(list)

    dataframe = np.array(array)

    ## Plotting data
    fig, ax = plt.subplots()
    im = ax.imshow(dataframe, cmap=color_palette)
    ax.set_title(f"% of proteins found in clusters (gap = {gap})")

    # Show all ticks and label them with the respective list entries
    ax.set_xticks(np.arange(len(labels)), labels=labels)
    ax.set_yticks(np.arange(len(labels)), labels=labels)

    # Create colorbar
    cbar_kw = {}
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel('', rotation=-90, va="bottom")

    # Tick labels
    plt.setp(
        ax.get_xticklabels(),
        rotation=45,
        ha="right",
        rotation_mode="anchor"
    )

    # Loop over data dimensions and create text annotations.
    for x in range(len(labels)):
        for y in range(len(labels)):
            text = ax.text(
                y,
                x,
                dataframe[x, y],
                ha="center",
                va="center",
                color=num_color
            )

    fig.tight_layout()

    filename = outdir + '/' + 'proteins_in_clusters.gap_' + gap + '.png'
    plt.savefig(filename)