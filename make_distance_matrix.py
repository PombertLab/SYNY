#!/usr/bin/python
## Pombert Lab 2022

name = 'make_distance_matrix.pl'
version = '0.3'
updated = '2022-06-17'

from os import listdir, mkdir
from os.path import isfile, basename, isdir
from sys import argv
import argparse
import re

usage = f"""
NAME		{name}
VERSION		{version}
UPDATED		{updated}
SYNOPSIS	Produces a distance matrix from syntetic clusters, pairs, and bidirectional
	homologs.

USAGE		{name} \\
		 -l SYNY_OUTPUT/LISTS \\
		 -c SYNY_OUTPUT/CLUSTERS

OPTIONS
-l (--lists)		Directory containing .list files
-n (--neighbors)	Directory containing .neighbors files
-p (--pairs)		Directory containing .pairs files
-c (--clusters)		Directory containing .clusters files
-o (--out)		Output directory
"""

def die(statement):
	print(statement)
	exit()

if len(argv) < 2:
	die(usage)

parser = argparse.ArgumentParser()
parser.add_argument("-l","--lists")
parser.add_argument("-n","--neighbors")
parser.add_argument("-p","--pairs")
parser.add_argument("-c","--clusters")
parser.add_argument("-o","--out",default="DIST_MAT")
args = parser.parse_args()

lists = args.lists
cluster_dir = args.clusters
pair_dir = args.pairs
neighbor_dir = args.neighbors
out = args.out

if not isdir(out):
	mkdir(out)

## Parse through list files getting the number of total proteins
names = []
genome_meta = {}
for file in listdir(lists):
	if isfile(f"{lists}/{file}"):
		with open(f"{lists}/{file}","r") as FILE:
			file_name = (basename(file)).replace(".list","")
			previous_accession = ""
			accessions = 0
			proteins = 0
			for line in FILE:
				line = line.replace("\n","")
				data = line.split("\t")
				accession = data[1]
				proteins += 1
				if accession != previous_accession:
					accessions += 1
				previous_accession = accession
			genome_meta[file_name] = [accessions,proteins]
			# print(f"{file_name[0:6]}\t{accessions}\t{proteins}")
			names.append(file_name)
names = sorted(names)

########################################################################################################################
## Generating a distance matrix from bidirectional searches
########################################################################################################################

neighbor_matrix = {}

for name_1 in names:
	for name_2 in names:
		if name_1 != name_2:

			neighbors = []

			protein_count = 0

			FILE = open(f"{neighbor_dir}/{name_1}_vs_{name_2}.neighbors","r")

			for line in FILE:
				line = line.replace("\n","")
				data = line.split("\t")
				neighbor_1 = data[0]
				neighbor_2 = data[2]\

				if neighbor_1 not in neighbors:
					neighbors.append(neighbor_1)
					protein_count += 1

				if neighbor_2 not in neighbors:
					neighbors.append(neighbor_2)
					protein_count += 1

			if name_1 not in neighbor_matrix:
				neighbor_matrix[name_1] = {}

			neighbor_matrix[name_1][name_2] = 1 - (protein_count/genome_meta[name_1][1])

		else:

			if name_1 not in neighbor_matrix:
				neighbor_matrix[name_1] = {}

			neighbor_matrix[name_1][name_2] = 0.0

for x,name_1 in enumerate(names):
	for y,name_2 in enumerate(names):
		if y <= x:
			value = (neighbor_matrix[name_1][name_2] + neighbor_matrix[name_2][name_1]) / 2
			neighbor_matrix[name_1][name_2] = value
			neighbor_matrix[name_2][name_1] = value

OUT = open(f"{out}/bidirectional_distance_matrix.tsv","w")
OUT.write(f"{names[0]}")
for val in names[1:]:
	OUT.write(f"\t{val}")
OUT.write("\n")
for name_1 in sorted(neighbor_matrix.keys()):
	key_set = [i for i in sorted(neighbor_matrix[name_1].keys())]
	OUT.write(f"{neighbor_matrix[name_1][key_set[0]]}")
	for name_2 in key_set[1:]:
		OUT.write(f"\t{neighbor_matrix[name_1][name_2]}")
	OUT.write("\n")
OUT.close()

########################################################################################################################
## Generating a distance matrix from pairs
########################################################################################################################

pairs_matrix = {}

for name_1 in names:
	for name_2 in names:
		if name_1 != name_2:

			neighbors = []

			pairs_count = 0

			FILE = open(f"{pair_dir}/{name_1}_vs_{name_2}.pairs","r")

			for line in FILE:
				line = line.replace("\n","")
				data = line.split("\t")
				neighbor_1 = data[0]
				neighbor_2 = data[2]\

				if neighbor_1 not in neighbors:
					neighbors.append(neighbor_1)
					pairs_count += 1

				if neighbor_2 not in neighbors:
					neighbors.append(neighbor_2)
					pairs_count += 1

			if name_1 not in pairs_matrix:
				pairs_matrix[name_1] = {}

			pairs_matrix[name_1][name_2] = 1 - (pairs_count/(genome_meta[name_1][1] - genome_meta[name_1][0]))

		else:

			if name_1 not in pairs_matrix:
				pairs_matrix[name_1] = {}

			pairs_matrix[name_1][name_2] = 0.0

for x,name_1 in enumerate(names):
	for y,name_2 in enumerate(names):
		if y <= x:
			value = (pairs_matrix[name_1][name_2] + pairs_matrix[name_2][name_1]) / 2
			pairs_matrix[name_1][name_2] = value
			pairs_matrix[name_2][name_1] = value

OUT = open(f"{out}/pair_distance_matrix.tsv","w")
OUT.write(f"{names[0]}")
for val in names[1:]:
	OUT.write(f"\t{val}")
OUT.write("\n")
for name_1 in sorted(pairs_matrix.keys()):
	key_set = [i for i in sorted(pairs_matrix[name_1].keys())]
	OUT.write(f"{pairs_matrix[name_1][key_set[0]]}")
	for name_2 in key_set[1:]:
		OUT.write(f"\t{pairs_matrix[name_1][name_2]}")
	OUT.write("\n")
OUT.close()

########################################################################################################################
## Generating a distance matrix from cluster files
########################################################################################################################

cluster_matrix = {}

for name_1 in names:
	for name_2 in names:
		if name_1 != name_2:
			protein_count = 0
			cluster_count = 0
			FILE = open(f"{cluster_dir}/{name_1}_vs_{name_2}.clusters","r")
			for line in FILE:
				line = line.replace("\n","")
				if "###" == line[0:3]:
					cluster_count += 1
				else:
					protein_count += 1

			if name_1 not in cluster_matrix:
				cluster_matrix[name_1] = {}

			if cluster_count == 0:
				cluster_matrix[name_1][name_2] = 1
			else:
				cluster_matrix[name_1][name_2] = 1 - ((protein_count/cluster_count)/(genome_meta[name_1][1]/genome_meta[name_1][0]))

		else:
			if name_1 not in cluster_matrix:
				cluster_matrix[name_1] = {}
			cluster_matrix[name_1][name_2] = 0.0

for x,name_1 in enumerate(names):
	for y,name_2 in enumerate(names):
		if y <= x:
			value = (cluster_matrix[name_1][name_2] + cluster_matrix[name_2][name_1]) / 2
			cluster_matrix[name_1][name_2] = value
			cluster_matrix[name_2][name_1] = value

OUT = open(f"{out}/cluster_distance_matrix.tsv","w")
OUT.write(f"{names[0]}")
for val in names[1:]:
	OUT.write(f"\t{val}")
OUT.write("\n")
for name_1 in sorted(cluster_matrix.keys()):
	key_set = [i for i in sorted(cluster_matrix[name_1].keys())]
	OUT.write(f"{cluster_matrix[name_1][key_set[0]]}")
	for name_2 in key_set[1:]:
		OUT.write(f"\t{cluster_matrix[name_1][name_2]}")
	OUT.write("\n")
OUT.close()

sum_of_all_matrix = {}

for name_1 in names:
	for name_2 in names:
		if name_1 not in sum_of_all_matrix.keys():
			sum_of_all_matrix[name_1] = {}
		sum_of_all_matrix[name_1][name_2] = neighbor_matrix[name_1][name_2]
		sum_of_all_matrix[name_1][name_2] += pairs_matrix[name_1][name_2]
		sum_of_all_matrix[name_1][name_2] += cluster_matrix[name_1][name_2]
		sum_of_all_matrix[name_1][name_2] /= 3

OUT = open(f"{out}/sum_of_all_distance_matrix.tsv","w")
OUT.write(f"{names[0]}")
for val in names[1:]:
	OUT.write(f"\t{val}")
OUT.write("\n")
for name_1 in sorted(sum_of_all_matrix.keys()):
	key_set = [i for i in sorted(sum_of_all_matrix[name_1].keys())]
	OUT.write(f"{sum_of_all_matrix[name_1][key_set[0]]}")
	for name_2 in key_set[1:]:
		OUT.write(f"\t{sum_of_all_matrix[name_1][name_2]}")
	OUT.write("\n")
OUT.close()