#!/usr/bin/env python3
## Pombert lab, 2024
version = '0.3a'
updated = '2024-05-27'
name = 'paf_metrics.py'

import os
import sys
import re
import argparse
import matplotlib.pyplot as plt
from multiprocessing import Pool, Value

################################################################################
## README
################################################################################

usage = f"""
NAME        {name}
VERSION     {version}
UPDATED     {updated}
SYNOPSIS    Calculates alignment metrics from PAF files and plots their 
            alignment lenghts vs. their similarity with matplotlib

COMMAND     {name} \\
              -f file.paf \\
              -c darkorange \\
              --threads 16

I/O OPTIONS:
-p (--paf)    PAF file to plot
-d (--outdir) Output directory [Default: ./]
-c (--color)  Color to use; red, green, blue... [Default: green]
              # https://matplotlib.org/stable/gallery/color/named_colors.html
-h (--height) Figure height in inches [Default: 10.8]
-w (--width)  Figure width in inches [Default: 19.2]
--fontsize    Figure font size [Default: 12]
--threads     Number of threads to use [Default: 16]
--version     Show script version
"""

# Print custom message if argv is empty
if (len(sys.argv) <= 1):
    print(usage)
    exit(0)

################################################################################
## Create command lines switches
################################################################################

cmd = argparse.ArgumentParser(add_help=False)
cmd.add_argument("-p", "--paf", nargs='*')
cmd.add_argument("-d", "--outdir", default='./')
cmd.add_argument("-c", "--color", default='steelblue')
cmd.add_argument("-h", "--height", default=10.8)
cmd.add_argument("-w", "--width", default=19.2)
cmd.add_argument("--fontsize", default=12)
cmd.add_argument("--threads", default=16)
cmd.add_argument("--version", action='store_true')
args = cmd.parse_args()

paf_files = args.paf
outdir = args.outdir
height = args.height
width = args.width
rgb = args.color
fontsize = int(args.fontsize)
threads = int(args.threads)
scversion = args.version

#########################################################################
### Version
#########################################################################

if scversion:
    print ("\b")
    print (f"Script:     {name}")
    print (f"Version:    {version}")
    print (f"Updated:    {updated}\n")
    exit(0)

################################################################################
## Function to calculate n metrics; e.g n50, n75, n90
################################################################################

def n_metric(list, n):
    n_threshold = int(sum(list)*n)
    n_sum = 0
    nmetric = 0

    for x in list:
        n_sum += x
        if n_sum >= n_threshold:
            nmetric = x
            break

    nmetric = "{:,}".format(nmetric)
    return nmetric

################################################################################
## Working on output directory
################################################################################

txtdir = outdir + '/TXT'
pngdir = outdir + '/PNG'
svgdir = outdir + '/SVG'

for dir in [outdir,txtdir,pngdir,svgdir]:
    if os.path.isdir(dir) == False:
        try:
            os.makedirs(dir)
        except:
            sys.exit(f"Can't create directory {dir}...")

################################################################################
## Working on PAF file
################################################################################

## Checking of PAF files are empty
non_zero_paf = []

for paf in paf_files:

    file_size = os.path.getsize(paf)
    basename = os.path.basename(paf)
    metrics_file = basename
    metrics_output = txtdir + '/' + metrics_file  + '.txt'

    if file_size == 0:
        warning = f"[WARNING]: paf_metrics.py - {paf} size is 0 byte."
        explanation = f"[WARNING]: paf_metrics.py - {paf} might have been terminated by OOM killer."
        print(warning, explanation, sep="\n", file=sys.stderr)

        metrics_output = txtdir + '/' + metrics_file
        txtout = open(metrics_output,'w')
        txtout.write(f"Metrics for {paf}:\n\n")
        txtout.write(f"PAF file is empty\n")

    else:
        non_zero_paf.append(paf)

## Plot if not empty
lsize = len(non_zero_paf)
counter = Value('i', 0)

def scatterplot(paf):

    global counter
    basename = os.path.basename(paf)
    metrics_file = basename

    m = re.search(r'(\S+)\.paf', basename)
    if m:
        prefix = m.group(1)
        metrics_file = prefix

    x_aln_sizes = []
    similarity = []
    matches = []

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
            n_matches = int(data[9])
            block_len = int(data[10])

            matches.append(n_matches)
            x_aln_sizes.append(block_len)

            sim = (n_matches/block_len)*100
            similarity.append(sim)

    ################################################################################
    ## Calculate alignment metrics
    ################################################################################

    aln_num = "{:,}".format(len(x_aln_sizes))
    aln_sum = "{:,}".format(sum(x_aln_sizes))
    longest_aln = "{:,}".format(max(x_aln_sizes))
    shortest_aln = "{:,}".format(min(x_aln_sizes))

    average = sum(x_aln_sizes)/len(x_aln_sizes)
    average = int(round(average))
    average = "{:,}".format(average)

    average_id = "{:.2f}".format( (sum(matches) / sum(x_aln_sizes)) * 100 )

    # n50, 75, 90
    ndata = x_aln_sizes[:]
    ndata.sort(reverse=True)
    n50 = n_metric(ndata,0.5)
    n75 = n_metric(ndata,0.75)
    n90 = n_metric(ndata,0.9)

    # median
    median_location = int(round(len(ndata)/2))
    median = ndata[median_location]
    median = "{:,}".format(median)

    # Print metrics to file or STDOUT
    metrics_output = txtdir + '/' + metrics_file  + '.txt'
    txtout = open(metrics_output,'w')

    txtout.write(f"Metrics for {paf}:\n\n")
    txtout.write(f"Bases + gaps in aligned blocks:\t{aln_sum}\n")
    txtout.write(f"Average sequence identity (%):\t{average_id}\n")
    txtout.write(f"# Alignments:\t{aln_num}\n")
    txtout.write(f"Longest:\t{longest_aln}\n")
    txtout.write(f"Shortest:\t{shortest_aln}\n")
    txtout.write(f"Average:\t{average}\n")
    txtout.write(f"Median:\t\t{median}\n")
    txtout.write(f"N50:\t\t{n50}\n")
    txtout.write(f"N75:\t\t{n75}\n")
    txtout.write(f"N90:\t\t{n90}\n")

    ################################################################################
    ## Plot alignment length vs sequence identity with matplotlib
    ################################################################################

    ##### Metrics text box
    adjust_l = 0
    metrics_list = [aln_sum, aln_num, longest_aln, shortest_aln, average, average_id, median, n50, n75, n90]
    for metric in metrics_list:
        if len(metric) > adjust_l:
            adjust_l = len(metric)

    metrics = f"""PAF metrics

        Bases + gaps in aligned blocks: {aln_sum.rjust(adjust_l + 1)}
        Average sequence identity (%): {average_id.rjust(adjust_l + 1)}
        # alignments: {aln_num.rjust(adjust_l + 1)}
        Longest: {longest_aln.rjust(adjust_l + 1)}
        Shortest: {shortest_aln.rjust(adjust_l + 1)}
        Average: {average.rjust(adjust_l + 1)}
        Median: {median.rjust(adjust_l + 1)}
        N50: {n50.rjust(adjust_l + 1)}
        N75: {n75.rjust(adjust_l + 1)}
        N90: {n90.rjust(adjust_l + 1)}"""

    metrics = metrics.replace("\t\t","")

    ##### Plotting data

    # Setting default image to widescreen by default
    plt.rcParams["figure.figsize"] = (width,height)
    plt.rcParams['svg.fonttype'] = 'none'
    plt.rcParams.update({'font.size': fontsize})

    # Setting labels
    xlabel = 'Alignment length (in Kb)'
    ylabel = 'Sequence identity (%)'

    title = os.path.basename(paf)
    title_font = 'normal'

    plt.title(
        title,
        loc='center',
        fontsize = fontsize,
        y = 0.2,
        pad=-50,
        fontweight=title_font
    )

    xmax = max(x_aln_sizes) / 1000
    plt.text(
        xmax,
        50,
        metrics,
        fontsize=fontsize,
        family='monospace',
        va='top',
        ha='right'
    )

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.ylim(0,105)
    plt.yscale('linear')
    plt.xscale('linear')

    plt.scatter(
        [x / 1000 for x in x_aln_sizes],
        similarity,
        color=rgb,
        alpha=0.75
    )

    # Creating outfiles
    with counter.get_lock():
        counter.value += 1
    pngfile = pngdir + '/' + metrics_file + '.png' 
    print(f"{counter.value} / {lsize*2} - plotting {pngfile}...")
    plt.savefig(pngfile)

    with counter.get_lock():
        counter.value += 1
    svgfile = svgdir + '/' + metrics_file + '.svg' 
    print(f"{counter.value} / {lsize*2} - plotting {svgfile}...")
    plt.savefig(svgfile)

    ## Close fig
    plt.clf()
    plt.cla()
    plt.close('all')

## Run
pool = Pool(threads)
pool.map(scatterplot, non_zero_paf)