#!/usr/bin/env python3
## Pombert lab, 2024

name = 'protein_cluster_hm.py'
version = '0.2c'
updated = '2024-03-30'

import sys
import os
import matplotlib
import matplotlib.pyplot as plt
import argparse
import re
import seaborn as sns
import pandas as pd

################################################################################
## README
################################################################################

usage = f"""
NAME        {name}
VERSION     {version}
UPDATED     {updated}
SYNOPSIS    Create heatmaps from SYNY matrices with matplotlib

REQS        matplotlib, pandas, seaborn

COMMAND    {name} \\
            -t clusters_summary_table.tsv \\
            -o HEATMAPS 

OPTIONS:
-t (--tsv)      Matrix in TSV format, e.g. matrix_gap_0.tsv
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

pngdir = outdir + '/PNG'
svgdir = outdir + '/SVG'

for dir in [outdir,pngdir,svgdir]:
    if os.path.isdir(dir) == False:
        try:
            os.makedirs(dir)
        except:
            sys.exit(f"Can't create directory {dir}...")

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

    clustered_png = pngdir + '/' + 'proteins_in_clusters.gap_' + gap + '.clustered.png'
    clustered_svg = svgdir + '/' + 'proteins_in_clusters.gap_' + gap + '.clustered.svg'

    heatmap_png = pngdir + '/' + 'proteins_in_clusters.gap_' + gap + '.heatmap.png'
    heatmap_svg = svgdir + '/' + 'proteins_in_clusters.gap_' + gap + '.heatmap.svg'

    ## Clustered heatmaps
    cm = sns.clustermap(
        data[0:],
        cmap=color_palette,
        annot=True,
        fmt='.1f'
    )

    cm.fig.suptitle(f"% of proteins found in clusters (gap = 0)", x=0.5, y=0.95)
    plt.savefig(clustered_png)
    plt.savefig(clustered_svg)
    plt.close()

    ## Normal heatmaps
    hm = sns.heatmap(
        data[0:],
        cmap=color_palette,
        annot=True,
        fmt='.1f'
    )

    hm.figure.suptitle(f"% of proteins found in clusters (gap = 0)", x=0.5, y=0.95)
    plt.savefig(heatmap_png)
    plt.savefig(heatmap_svg)
    plt.close()