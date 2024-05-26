#!/usr/bin/env python3
## Pombert lab, 2024

version = '0.2b'
updated = '2024-05-12'
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
SYNOPSIS    Checks available matplotlib + seaborn color palettes on the system

REQS        matplotlib, seaborn

COMMAND    {name} \\
             --list \\
             --plot \\
             --outfile palettes.png
OPTIONS:
-l (--list)     List available color palettes
-p (--plot)     Plot available color palettes
-o (--outfile)  Output file [Default: color_palettes.png]
-h (--height)   Figure height in inches [Default: 10.8]
-w (--width)    Figure width in inches [Default: 19.2]
-f (--font)     Font size [Default: 6]
-c (--check)    Check if color palettes entered exists
"""

# Print custom message if argv is empty
if (len(sys.argv) <= 1):
    print(usage)
    exit(0)

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

if colorlist:

    print(f"\n##### Available color palettes #####\n")

    ## matplotlib
    mcolors = list(sorted(colormaps))
    numcol = len(mcolors)

    ## matplotlib + seaborn
    import seaborn as sns
    matsea = list(sorted(colormaps))
    catcol = len(matsea)

    ## seaborn only
    seaborn = []
    for x in matsea:
        if x not in mcolors:
            seaborn.append(x)
    seacol = len(seaborn)

    ## print color info
    print(f"# matplotlib + seaborn ({catcol}):\n")
    print(', '.join(matsea), "\n")

    print(f"# matplotlib only ({numcol}):\n")
    print(', '.join(mcolors), "\n")

    print(f"# seaborn only ({seacol}):\n")
    print(', '.join(seaborn), "\n")

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
    
    if len(false_colors) == 0:
        print("\n[OK] All requested color palettes are available:", *check, "\n")
        exit()

################################################################################
## Plot available color palettes
################################################################################

import seaborn as sns

mcolors = list(sorted(colormaps))

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