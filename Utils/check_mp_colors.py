#!/usr/bin/env python3
## Pombert lab, 2024

version = '0.2a'
updated = '2024-05-07'
name = 'check_mp_colors.py'

import argparse
from matplotlib import colormaps
from pylab import *
from numpy import outer

################################################################################
## README
################################################################################

usage = f"""
NAME        {name}
VERSION     {version}
UPDATED     {updated}
SYNOPSIS    Checks available matplotlib + seaborn colors on the system

REQS        matplotlib, seaborn

COMMAND    {name} \\
             --list \\
             --plot \\
             --outfile palettes.png
OPTIONS:
-l (--list)     Lists available colors
-p (--plot)     Plots available colors
-o (--outfile)  Output file [Default: color_palettes.png]
-h (--height)   Figure height in inches [Default: 10.8]
-w (--width)    Figure width in inches [Default: 19.2]
-f (--font)     Font size [Default: 6]
-c (--check)    Check if colors entered exists
"""

# Print custom message if argv is empty
if (len(sys.argv) <= 1):
    print(usage)
    sys.exit()

################################################################################
## Create command lines switches
################################################################################

cmd = argparse.ArgumentParser(add_help=False)
cmd.add_argument("-l", "--list", action='store_true')
cmd.add_argument("-p", "--plot", action='store_true')
cmd.add_argument("-o", "--outfile", default='color_palettes.png')
cmd.add_argument("-w", "--width", default=60)
cmd.add_argument("-h", "--height", default=6)
cmd.add_argument("-f", "--font", default=12)
cmd.add_argument("-c", "--check", nargs='*')
args = cmd.parse_args()

colorlist = args.list
colorplot = args.plot
outfile = args.outfile
height = args.height
width = args.width
fontsize = args.font
check = args.check

################################################################################
## Check for available color palettes
################################################################################

mcolors = list(sorted(colormaps))

if colorlist:

    print(f"\n### Available color palettes\n")

    ## matplotlib
    mcolors = list(sorted(colormaps))
    numcol = len(mcolors)
    print(f"matplotlib only ({numcol}):\n")
    print(', '.join(mcolors), "\n")

    ## matplotlib + seaborn
    import seaborn as sns
    mcolors = list(sorted(colormaps))
    numcol = len(mcolors)
    print(f"matplotlib + seaborn ({numcol}):\n")
    print(', '.join(mcolors), "\n")

################################################################################
## Perform check for colors then exit
################################################################################

if check:

    import seaborn as sns
    mcolors = list(sorted(colormaps))

    false_colors = []
    dcolors = {}

    for palette in mcolors:
        dcolors[palette] = 1

    for palette in check:
        if palette in dcolors:
            next
        else:
            false_colors.append(palette)

    if len(false_colors) >= 1:
        print("\n[E] The following color palettes are not available:", *false_colors)
        print(f"\nAvailable matplotlib + seaborn color palettes ({len(mcolors)}):\n")
        print(', '.join(mcolors), "\n")
        exit()

################################################################################
## Plot available color palettes
################################################################################

import seaborn as sns

plt.rcParams["figure.figsize"] = (width,height)

if colorplot:

    a = outer(arange(0,1,0.01),ones(10))

    subplots_adjust(
        top=0.9,
        bottom=0.05,
        left=0.01,
        right=0.99
    )

    palette_num = len(mcolors) + 1

    for x, name in enumerate(mcolors):
        subplot(1, palette_num, x + 1)
        axis("off")
        imshow(a, aspect='equal', cmap=get_cmap(name), origin='lower')
        title(name, rotation=90, fontsize=fontsize)

    plt.savefig(outfile)