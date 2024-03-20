#!/usr/bin/env python3
## Pombert lab, 2024
version = '0.1'
updated = '2024-03-20'
name = 'paf_metrics.py'

import os
import sys
import argparse
import matplotlib.pyplot as plt

################################################################################
## README
################################################################################

usage = f"""
NAME		{name}
VERSION		{version}
UPDATED		{updated}
SYNOPSIS	Calculates alignment metrics from PAF files and plots their 
			alignment lenghts vs. their similarity with matplotlib

COMMAND		{name} \\
		  -f file.paf \\
		  -c darkorange \\
		  -o aln_distribution.svg aln_distribution.pdf

I/O OPTIONS:
-p (--paf)	PAF file to plot
-d (--outdir)	Output directory [Default: ./]
-m (--metrics)	Metrics output file [Default: aln_metrics.txt]
-o (--output)	Save plot to specified output file(s) (supported formats: png, pdf, svg)

PLOT OPTIONS:
-c (--color)	Color to use; red, green, blue... [Default: green]
		# https://matplotlib.org/stable/gallery/color/named_colors.html
-h (--height)	Figure height in inches [Default: 10.8]
-w (--width)	Figure width in inches [Default: 19.2]
"""

# Print custom message if argv is empty
if (len(sys.argv) <= 1):
	print(usage)
	sys.exit()

################################################################################
## Create command lines switches
################################################################################

cmd = argparse.ArgumentParser(add_help=False)
cmd.add_argument("-p", "--paf")
cmd.add_argument("-o", "--output", nargs='*')
cmd.add_argument("-d", "--outdir", default='./')
cmd.add_argument("-m", "--metrics", default='aln_metrics.txt')
cmd.add_argument("-c", "--color", default='steelblue')
cmd.add_argument("-h", "--height", default=10.8)
cmd.add_argument("-w", "--width", default=19.2)
args = cmd.parse_args()

paf = args.paf
output = args.output
outdir = args.outdir
metrics_file = args.metrics
height = args.height
width = args.width
rgb = args.color

################################################################################
## Working on output directory
################################################################################

if output is not None:
	if os.path.isdir(outdir) == False:
		try:
			os.makedirs(outdir)
		except:
			sys.exit(f"Can't create directory {outdir}...")


################################################################################
## Working on PAF file
################################################################################

x_aln_sizes = []
similarity = []

FH = open(paf,'r')
print(f"\nCalculating alignment metrics for: {paf} ...")

## Parse alignments
for line in FH:

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

    x_aln_size = block_len
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

# Function to calculate n metrics; e.g n50, n75, n90
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
pmetrics = f"""Metrics for {paf}:

Bases + gaps in aligned blocks:\t{aln_sum}
# Alignments:\t{aln_num}
Longest:\t{longest_aln}
Shortest:\t{shortest_aln}
Average:\t{average}
Median:\t\t{median}
N50:\t\t{n50}
N75:\t\t{n75}
N90:\t\t{n90}
"""

metrics_output = outdir + '/' + metrics_file
METRICS = open(metrics_output,'w')
print(pmetrics, file=METRICS)

################################################################################
## Plot alignment length vs sequence identity with matplotlib
################################################################################

##### Metrics text box

adjust_l = 0
metrics_list = [aln_sum, aln_num, longest_aln, shortest_aln, average, median, n50, n75, n90]
for metric in metrics_list:
	if len(metric) > adjust_l:
		adjust_l = len(metric)

metrics = f"""
	PAF metrics

	Bases + gaps in aligned blocks: {aln_sum.rjust(adjust_l + 1)}
	# alignments: {aln_num.rjust(adjust_l + 1)}
	Longest: {longest_aln.rjust(adjust_l + 1)}
	Shortest: {shortest_aln.rjust(adjust_l + 1)}
	Average: {average.rjust(adjust_l + 1)}
	Median: {median.rjust(adjust_l + 1)}
	N50: {n50.rjust(adjust_l + 1)}
	N75: {n75.rjust(adjust_l + 1)}
	N90: {n90.rjust(adjust_l + 1)}
"""
metrics = metrics.replace("\t","")

##### Plotting data

# Setting default image to widescreen by default
plt.rcParams["figure.figsize"] = (width,height)

# Setting labels
xlabel = 'Alignment length (in Kb)'
ylabel = 'Sequence identity (%)'

title = os.path.basename(paf)
title_font = 'normal'

plt.title(
	title,
	loc='center',
	fontsize = 12,
	y = 0.2,
	pad=-50,
	fontweight=title_font
)

xmax = max(x_aln_sizes) / 1000
ymax = max(similarity)

plt.text(
	xmax,
	50,
	metrics,
	fontsize=10,
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

# Creating outfile
for x in output:
    filename = outdir + '/' + x 
    print(f"Plotting {filename}...")
    plt.savefig(filename)
