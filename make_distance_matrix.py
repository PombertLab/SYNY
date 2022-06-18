#!/usr/bin/python
## Pombert Lab 2022

name = 'make_distance_matrix.pl'
version = '0.2'
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
-c (--clusters)		Directory containing .clusters files
-p (--pairs)		Directory containing .pairs files
-o (--out)		Output directory
"""

def die(statement):
	print(statement)
	exit()

if len(argv) < 2:
	die(usage)

parser = argparse.ArgumentParser()
parser.add_argument("-l","--lists")
parser.add_argument("-c","--clusters")
parser.add_argument("-p","--pairs")
parser.add_argument("-o","--out",default="DIST_MAT")
args = parser.parse_args()

lists = args.lists
clusters = args.clusters
pairs = args.pairs
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
			print(f"{file_name[0:6]}\t{accessions}\t{proteins}")
			names.append(file_name)
names = sorted(names)

cluster_matrix = {}

for name_1 in names:
	for name_2 in names:
		if name_1 != name_2:
			protein_count = 0
			cluster_count = 0
			FILE = open(f"{clusters}/{name_1}_vs_{name_2}.clusters","r")
			for line in FILE:
				line = line.replace("\n","")
				if "###" == line[0:3]:
					cluster_count += 1
				else:
					protein_count += 1
			if name_1 not in cluster_matrix:
				cluster_matrix[name_1] = {}
			cluster_matrix[name_1][name_2] = 1 - (protein_count/genome_meta[name_1][1])
			# if cluster_count != 0:
			# 	cluster_matrix[name_1][name_2] = 1 - ((protein_count/cluster_count)/(genome_meta[name_1][1]/genome_meta[name_1][0]))
			# else:
			# 	cluster_matrix[name_1][name_2] = 1
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

for x,name_1 in enumerate(names):
	print(f"{name_1[0:6]}",end="\t")
	for y,name_2 in enumerate(names):
		print(round(cluster_matrix[name_1][name_2],3),end="\t")
	print()

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