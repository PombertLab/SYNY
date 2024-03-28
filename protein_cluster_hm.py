#!/usr/bin/env python3
## Pombert lab, 2024

name = 'protein_cluster_hm.py'
version = '0.2a'
updated = '2024-03-28'

import sys
import os
import matplotlib
import matplotlib.pyplot as plt
import argparse
import re
import numpy as np
import seaborn as sns
import pandas as pd

################################################################################
## README
################################################################################

usage = f"""
NAME        {name}
VERSION     {version}
UPDATED     {updated}
SYNOPSIS    Create heatmaps from SYNY summary files with matplotlib

REQS        matplotlib, pandas, seaborn

COMMAND    {name} \\
            -t clusters_summary_table.tsv \\
            -o HEATMAPS 

OPTIONS:
-t (--tsv)      clusters_summary_table.tsv
-o (--outdir)   Output directory [Default: ./HEATMAPS]
-p (--palette)  Color palette [Default: crest]
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
args = cmd.parse_args()

tsv_file = args.tsv
outdir = args.outdir
height = args.height
width = args.width
color_palette = args.palette

################################################################################
## Working on output directory
################################################################################

if os.path.isdir(outdir) == False:
    try:
        os.makedirs(outdir)
    except:
        sys.exit(f"Can't create directory {outdir}...")

################################################################################
## ### 
################################################################################

plt.rcParams["figure.figsize"] = (width,height)

with open (tsv_file) as f:

    data = pd.read_csv(tsv_file, sep="\t", header=0, index_col=0)

    gap = None
    m = re.search(r'gap_(\d+).tsv', tsv_file)
    if m:
        gap = m.group(1)

    clustered = outdir + '/' + 'proteins_in_clusters.gap_' + gap + '.clustered'
    heatmap = outdir + '/' + 'proteins_in_clusters.gap_' + gap + '.heatmap'

    cm = sns.clustermap(
        data[0:],
        cmap=color_palette,
        annot=True,
        fmt='.1f'
    )

    cm.figure.suptitle(f"% of proteins found in clusters (gap = 0)", x=0.5, y=0.95)
    plt.savefig(f"{clustered}.png")
    plt.savefig(f"{clustered}.svg")
    plt.close()

    hm = sns.heatmap(
        data[0:],
        cmap=color_palette,
        annot=True,
        fmt='.1f'
    )

    hm.figure.suptitle(f"% of proteins found in clusters (gap = 0)", x=0.5, y=0.95)
    plt.savefig(f"{heatmap}.png")
    plt.savefig(f"{heatmap}.svg")
    plt.close()