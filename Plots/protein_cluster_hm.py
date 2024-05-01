#!/usr/bin/env python3
## Pombert lab, 2024

name = 'protein_cluster_hm.py'
version = '0.3c'
updated = '2024-05-01'

import sys
import os
import matplotlib.pyplot as plt
import argparse
import re
import seaborn as sns
import pandas as pd
from multiprocessing import Pool, Value

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
-p (--palette)  Color palette [Default: winter_r]
-h (--height)   Figure height in inches [Default: 10]
-w (--width)    Figure width in inches [Default: 10]
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
cmd.add_argument("-t", "--tsv", nargs='*')
cmd.add_argument("-o", "--outdir", default='./HEATMAPS')
cmd.add_argument("-h", "--height", default=10)
cmd.add_argument("-w", "--width", default=10)
cmd.add_argument("-p", "--palette", default='winter_r')
cmd.add_argument("--fontsize", default=8)
cmd.add_argument("--threads", default=16)
args = cmd.parse_args()

tsv_files = args.tsv
outdir = args.outdir
height = args.height
width = args.width
color_palette = args.palette
fontsize = int(args.fontsize)
threads = int(args.threads)

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
## Plotting heatmaps 
################################################################################

lsize = len(tsv_files) * 4
counter = Value('i', 0)

def heatmap(tsv_file):

    plt.rcParams["figure.figsize"] = (width,height)
    plt.rcParams['svg.fonttype'] = 'none'
    plt.rcParams.update({'font.size': fontsize})
    global counter

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
            fmt='.1f',
            vmin=0,
            vmax=100,
        )

        cm.fig.suptitle(f"% of proteins found in clusters (gap = {gap})", x=0.5, y=0.95)
        with counter.get_lock():
            counter.value += 1
        print(f"{counter.value} / {lsize} - plotting {clustered_png}")
        plt.savefig(clustered_png)

        with counter.get_lock():
            counter.value += 1
        print(f"{counter.value} / {lsize} - plotting {clustered_svg}")
        plt.savefig(clustered_svg)
        plt.close('all')

        ## Normal heatmaps
        hm = sns.heatmap(
            data[0:],
            cmap=color_palette,
            annot=True,
            fmt='.1f',
            vmin=0,
            vmax=100
        )

        hm.figure.suptitle(f"% of proteins found in clusters (gap = {gap})", x=0.5, y=0.95)
        with counter.get_lock():
            counter.value += 1
        print(f"{counter.value} / {lsize} - plotting {heatmap_png}")
        plt.savefig(heatmap_png)

        with counter.get_lock():
            counter.value += 1
        print(f"{counter.value} / {lsize} - plotting {heatmap_svg}")
        plt.savefig(heatmap_svg)

        plt.close('all')

## Run
pool = Pool(threads)
pool.map(heatmap, tsv_files)