#!/usr/bin/env python3
## Pombert lab, 2024
version = '0.1b'
updated = '2024-05-01'
name = 'paf_to_heatmap.py'

import sys
import os
import re
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

################################################################################
## README
################################################################################

usage = f"""
NAME        {name}
VERSION     {version}
UPDATED     {updated}
SYNOPSIS    Create heatmaps from PAF files with matplotlib

REQS        matplotlib, seaborn

COMMAND    {name} \\
            -p *.paf \\
            -o HEATMAPS \\
            -x matrix.txt

OPTIONS:
-p (--paf)      PAF file(s) to plot
-f (--fasta)    FASTA files used to generate PAF alignments
-o (--outdir)   Output directory [Default: ./]
-h (--height)   Figure height in inches [Default: 10]
-w (--width)    Figure width in inches [Default: 10]
-c (--palette)  Seaborn color palette [Default: winter_r]
                # See https://www.practicalpythonfordatascience.com/ap_seaborn_palette
                # for a list of color palettes
--fontsize      Font size [Default: 8]
-x (--matrix)   Matrix output file [Default: matrix.tsv]
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
cmd.add_argument("-h", "--height", default=10)
cmd.add_argument("-w", "--width", default=10)
cmd.add_argument("-c", "--palette", default='winter_r')
cmd.add_argument("-x", "--matrix", default='matrix.tsv')
cmd.add_argument("--fontsize", default=8)
args = cmd.parse_args()

paf_files = args.paf
fasta_files = args.fasta
outdir = args.outdir
height = args.height
width = args.width
color_palette = args.palette
fontsize = int(args.fontsize)
matrix_file = args.matrix

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
        m = re.search(r'(\S+)\.fasta', basename)
        if m:
            basename = m.group(1)
            len_dict[basename] = 0

        with open(fasta) as f:
            for line in f:
                if '>' not in line:
                    length = len(line)
                    len_dict[basename] += length

### Matrix dataframe
matrix = {}
for query in len_dict:
    matrix[query] = {}
    for subject in len_dict:
        matrix[query][subject] = 0
        if query == subject:
            matrix[query][subject] = 100

################################################################################
## Parsing and plotting PAF file(s) 
################################################################################

for paf in paf_files:

    basename = os.path.basename(paf)
    qfile = None
    sfile = None

    m = re.search(r'^(\w+)_vs_(\w+)', basename)
    if m:
        sfile = m.group(1)
        qfile = m.group(2)

    with open(paf) as file:
    
        dataframe = {}

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

            if query not in dataframe:
                dataframe[query] = []

            pair = str(q_start) + ',' + str(q_end)
            dataframe[query].append(pair)

        msum = 0
        ## Working per contig (query) to reduce mem usage
        for query in dataframe:

            overlaps = []

            for pair in dataframe[query]:
                data = pair.split(',')
                for x in range(int(data[0]), int(data[1]) + 1):
                    overlaps.append(x)

            ## Using set to remove overlapping bases
            qsetlen = len(set(overlaps))
            msum += qsetlen

        percentage = (msum / len_dict[qfile]) * 100
        matrix[qfile][sfile] = percentage

with open(matrix_file, "w") as file:

    for query in matrix:
        file.write(f"\t{query}")
    file.write("\n")

    for query in matrix:
        file.write(query)
        for subject in matrix[query]:
            file.write(f"\t{'{:0.2f}'.format(matrix[query][subject])}")
        file.write("\n")

################################################################################
## Plotting heatmaps
################################################################################

plt.rcParams["figure.figsize"] = (width,height)
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams.update({'font.size': fontsize})

with open (matrix_file) as f:

    data = pd.read_csv(matrix_file, sep="\t", header=0, index_col=0)

    clustered_png = pngdir + '/' + 'colinear_bases.mmap.clustered.png'
    clustered_svg = svgdir + '/' + 'colinear_bases.mmap.clustered.svg'

    heatmap_png = pngdir + '/' + 'colinear_bases.mmap.heatmap.png'
    heatmap_svg = svgdir + '/' + 'colinear_bases.mmap.heatmap.svg'

    ## Clustered heatmaps
    cm = sns.clustermap(
        data[0:],
        cmap=color_palette,
        annot=True,
        fmt='.1f',
        vmin=0,
        vmax=100
    )

    cm.fig.suptitle(f"% of colinear bases in pairwise alignments", x=0.5, y=0.95)
    print(f"1 / 4 - Plotting {clustered_png}")
    plt.savefig(clustered_png)

    print(f"2 / 4 - Plotting {clustered_svg}")
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

    hm.figure.suptitle(f"% of colinear bases in pairwise alignments", x=0.5, y=0.95)
    print(f"3 / 4 - Plotting {heatmap_png}")
    plt.savefig(heatmap_png)

    print(f"4 / 4 - Plotting {heatmap_svg}")
    plt.savefig(heatmap_svg)

    plt.close('all')