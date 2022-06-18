#!/usr/bin/python
## Pombert Lab 2022

name = 'phylogenetic_figure_maker.py'
version = '0.2'
updated = '2022-06-17'

import argparse as ap
from sys import argv
from os import mkdir
from os.path import isdir
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec

usage = f"""
NAME		{name}
VERSION		{version}
UPDATED		{updated}
SYNOPSIS	Uses neighbor-joining to produce a phylogentic tree from a distance matrix, 
		and creates a figure containing said tree in tandem with the distance matrix.

USAGE		{name} \\
		  -d SYNY/SYNTENY/gap_0/distance_matrix.tsv

OPTIONS
-d (--dist_mat)		Distance matrix (see GUIDE for file configuration)
-p (--prefix)		Output file prefix [Default: PHYLOGENTIC_FIGURE]
-f (--format)		Figure format(s) (.png,.eps,.jpg,.jpeg,.pdf,.pgf,.ps,.raw,.rgba,.svg,
			.svgz,.tif,.tiff) [Default: .svg]
-o (--out)		Directory containing figures [Default: PHYLOGENTIC_FIGURES]

GUIDE

The input file is a tab-delimited.
The first line should contain the organisms names.
Each following line should contain a row of the distance matrix.

Example:

E_coli	L_mono	S_cere	S_pombe	S_suis
0.0	0.952	1.0	1.0	0.961
0.952	0.0	1.0	0.999	0.856
1.0	1.0	0.0	0.997	1.0
1.0	0.999	0.997	0.0	1.0
0.961	0.856	1.0	1.0	0.0
"""

def die(statement):
	print(statement)
	exit()

if len(argv) < 2:
	die(usage)

class Tree:

	## Setting up initial variables and calling neighbor-joining and matplotlib functions
	def __init__(self,matrix,names,prefix,formats,outdir):

		## Distance matrix used for neighbor-joining
		self.Dij = matrix
		## List of names
		self.names = names
		## Output file prefix
		self.prefix = prefix
		## Figure save formats
		self.formats = formats
		## Output directory
		self.outdir = outdir

		## This will be a list of the branches
		self.tree = ""
		## This will be the LUCA node of the tree
		self.root = ""
		## List containing the non-joined species of the tree
		self.leaves = []

		## Generate the leaves
		for index,name in enumerate(names):
			self.__create_leaf(name,index)
		
		## Compress the leaves into a main root via neighbor-joinging
		self.__find_root(self.Dij,self.leaves,0)

		self.__create_figure()

	## Leaf class (OG species that are being neighbor-joined)
	class Leaf:

		## Setting up initial varialbes
		def __init__(self,name,index):
			## Leaf name
			self.name = name
			## The location of the leaf in the original list
			self.index = index
			## The parent node that the leaf belongs to (will be set during neighbor-joining)
			self.parent = False
			## The leaf/node that the leaf is joined to (will also be set during neighbor-joining)
			self.sibling = False
			## The level of joining this leaf resides at (this is used soley for creating the phylogenetic tree)
			self.level = ""
		
		## Sets parent and sibiling during neighbor-joining
		def add_connections(self,parent,sibling):
			self.parent = parent
			self.sibling = sibling

		## Sets the joining level
		def add_level(self,level):
			self.level = level

	## Node class (joining points of leafs/nodes)
	class Node:

		## Initialize variables (similar to leaf class, but contains branch information)
		def __init__(self,name,branch1,length1,branch2,length2,level):
			## The name of the node
			self.name = name
			## The two leafs/nodes being joined in this node
			self.branches = [branch1,branch2]
			## The length of the leafs/nodes to the joining node, respectively
			self.branch_length = [length1,length2]
			## The parent node that the node belongs to (will be set during neighbor-joining)
			self.parent = False
			## The leaf/node that the node is joined to (will also be set during neighbor-joining)
			self.sibling = False
			## The level of joining this node resides at (this is used soley for creating the phylogenetic tree)
			self.level = level

		## Sets parent and sibiling during neighbor-joining
		def add_connections(self,parent,sibling):
			self.parent = parent
			self.sibling = sibling

	## Root of the tree. Contains the LUCA of the provided species
	class Root:

		## Initialize variables (similar to node class)
		def __init__(self,name,branch1,length1,branch2,length2,level):
			self.name = name
			self.branches = [branch1,branch2]
			self.branch_length = [length1,length2]
			self.level = level

	## Creates a branch, connecting two leafs
	def __create_branch(self,name,branch1,length1,branch2,length2,level):
		return self.Node(name,branch1,length1,branch2,length2,level)

	## Go through creating the figure
	def __create_figure(self):

		def add_colorbar_outside(im,ax):
			fig = ax.get_figure()
			bbox = ax.get_position() #bbox contains the [x0 (left), y0 (bottom), x1 (right), y1 (top)] of the axis.
			width = 0.01
			eps = 0.01 #margin between plot and colorbar
			# [left most position, bottom position, width, height] of color bar.
			cax = fig.add_axes([bbox.x1 + eps, bbox.y0, width, bbox.height])
			cbar = fig.colorbar(im, cax=cax)

		def graph_tree(tree,array,plot,offsets):

			for val in array:
				if isinstance(val,list):
					graph_tree(tree,val,plot,offsets)

			if not isinstance(array[0],list) and not isinstance(array[1],list):
				
				index_a,xoffset_a,yoffset_a = offsets[array[0].name]
				index_b,xoffset_b,yoffset_b = offsets[array[1].name]

				plot.plot([-xoffset_a,-array[0].parent.level],[-index_a,-index_a],solid_capstyle='round',color='k')
			
				plot.plot([-xoffset_b,-array[1].parent.level],[-index_b,-index_b],solid_capstyle='round',color="k")
				
				plot.plot([-array[0].parent.level,-array[0].parent.level],[-index_a,-index_b],solid_capstyle='round',color="k")

				offsets[array[0].parent.name] = [(index_a+index_b)/2,array[0].parent.level,(index_a+index_b)/2]

				self.__exchange_siblings_for_parent(array,tree,array[0].parent)

		## Get the order of names from the newly formed tree
		def order_names(array,order):
			
			for i,val in enumerate(array):
				if isinstance(val,list):
					order_names(val,order)
				else:
					order.append(val.name)
		
		## Get the order of the species names
		name_order = []
		order_names([i for i in self.tree],name_order)
		
		## Prepare the line offsets for the tree
		graph_offsets = {name:[i,0,0] for i,name in enumerate(name_order)}
		
		## Get a distance matrix that is reorder to match the order of the tree
		new_mat = self.__map_data(self.Dij,self.names,name_order)

		## Figure object
		fig = plt.figure()

		## GridSpec object helps with organizing the plots
		## 1 row by 2 columns
		## no white space between the two
		## The second column (matrix) being twice the aspect ratio of the first (tree)
		gs = GridSpec(1,2,figure=fig,wspace=0.0,width_ratios=[1,2])#,height_ratios=[1,1])

		## Assign the tree the first position
		ax = fig.add_subplot(gs[0])

		## Assign the matrix the second position
		bx = fig.add_subplot(gs[1])

		im = bx.imshow(new_mat,aspect='equal',cmap='magma')
		bx.set_yticks([i for i in range(len(graph_offsets))])
		bx.set_yticklabels(graph_offsets)
		bx.set_xticks([i for i in range(len(graph_offsets))])
		bx.set_xticklabels(graph_offsets,rotation=90)
		bx.xaxis.tick_top()

		graph_tree(self.tree,[self.tree],ax,graph_offsets)
		ax.axis('off')
		gs.tight_layout(fig,pad=5,w_pad=0)
		add_colorbar_outside(im,bx)

		for form in self.formats:
			plt.savefig(f"{self.outdir}/{self.prefix}.{form}",format=form,transparent=True)

		# plt.show()

	## Creates a leaf
	def __create_leaf(self,name,index):
		self.leaves.append(self.Leaf(name,index))

	## Decompresses the node into child nodes/leafs
	def __exchange_siblings_for_parent(self,value,array,new_val):
		for i,val in enumerate(array):
			if val == value:
				array[i] = new_val
			elif type(val) == list:
				self.__exchange_siblings_for_parent(value,val,new_val)

	## Neighbor-joining
	def __find_root(self,Dij,leaves,node_counter):

		def get_level(node1,node2):
			if isinstance(node1,Tree.Leaf) and isinstance(node2,Tree.Leaf):
				level = 1
				node1.add_level(level)
				node2.add_level(level)
			elif isinstance(node1,Tree.Leaf):
				level = node2.level + 1
				node1.add_level(level)
			elif isinstance(node2,Tree.Leaf):
				level = node1.level + 1
				node2.add_level(level)
			else:
				level = max(node1.level,node2.level) + 1
			return level

		## Additive row distance
		S = [sum(Dij[i]) for i in range(len(Dij))]

		## Modified distance matrix
		Mij = [[0 for i in range(len(Dij[j]))] for j in range(len(Dij))]

		## The minimum value in the distance matrix
		min_val = 0
		## The indicies that the minimum value is located
		x,y = 0,0
		## The bramch distances of the minimum value
		xl,yl = 0,0
		
		## Iterate over the distance matrix to determine the minimum value and fill the modified distance matrix
		for i in range(len(Dij)):
			for j in range(len(Dij[i])):

				## The distance between a leaf/node and itself will always be the minimum value, so just leave it at 0
				if i != j:
					
					## Calculate the modified distance matrix value
					Mij[i][j] = ((len(Dij)-2)*Dij[i][j]) - (S[i] + S[j])
					
					## If the current modified distance value is the lowest
					if Mij[i][j] < min_val:
						
						min_val = Mij[i][j]
						x,y = i,j

						## Calculate the branch lengths between the two leaf/nodes and the joining node
						xl = (Mij[i][j]*.5) + ((1/2)*(len(Dij)-2)*abs(S[i] - S[j]))
						yl = Mij[i][j] + xl

		## New distance matrix that takes into account the joining of two leaf/nodes
		NDij = [[0 for i in range(len(Dij[j])-1)] for j in range(len(Dij)-1)]

		## Iterate over Dij, skipping the index belonging to the joined values
		x_offset = 0
		for i in range(len(Dij)):
			y_offset = 0
			if i == x or i == y:
				x_offset += 1
			else:
				for j in range(len(Dij[i])):
					if j == x or j == y:
						y_offset += 1
					else:
						NDij[i-x_offset+1][j-y_offset+1] = Dij[i][j]
		
		## Iterate over Dij, only filling index belonging to joined values
		x_offset = 0
		for i in range(len(Dij)):
			if i == x or i == y:
				x_offset += 1
			elif i < len(Dij):
				NDij[i-x_offset+1][0] = (Dij[x][i]+Dij[y][i]-Dij[x][y])/2
				NDij[0][i-x_offset+1] = (Dij[x][i]+Dij[y][i]-Dij[x][y])/2

		## Modified leaves to put back into __find_root()
		new_leaves = []

		## Create a new node that contains the two joined leaf/nodes
		new_node = self.__create_branch(f"U{node_counter}",leaves[x],xl,leaves[y],yl,get_level(leaves[x],leaves[y]))
		
		## Assign parent/sibling values to the joined leafs/nodes
		leaves[x].add_connections(new_node,leaves[y])
		leaves[y].add_connections(new_node,leaves[x])

		## Add the new node to the modified leaves list
		new_leaves.append(new_node)

		## Add the remaining leaves/nodes to the modified leaves list
		for index,leaf in enumerate(leaves):
			if index != x and index != y:
				new_leaves.append(leaf)

		## Increase the node number
		node_counter += 1

		## If the number of leaves/nodes in the modified distance matrix is greater 1, the root has not been determined,
		## rerun __find_root() with new values
		if len(NDij) > 1:
			tree = self.__find_root(NDij,new_leaves,node_counter)
			return tree
		## If the root as been reached, grow out the tree with the given relations
		else:
			self.root = new_node
			self.tree = new_node.branches
			self.__grow_tree(self.tree)

	## Takes tree root and expands out to leaves, producing the tree
	def __grow_tree(self,array):

		for val in array:
			if isinstance(val,Tree.Node):
				self.__exchange_siblings_for_parent(val,array,val.branches)
				self.__grow_tree(array)
			elif isinstance(val,list):
				self.__grow_tree(val)

	## Map original distance matrix to one that matches the order of the tree
	def __map_data(self,old_mat,old_names,new_names):
		
		new_mat = [[0 for i in range(len(Dij))] for j in range(len(Dij))]
		
		matrix = {}

		for x,row in enumerate(old_mat):
			x_name = old_names[x]
			if x_name not in matrix.keys():
				matrix[x_name] = {}
			for y,val in enumerate(row):
				y_name = old_names[y]
				matrix[x_name][y_name] = val

		for x,x_name in enumerate(new_names):
			for y,y_name in enumerate(new_names):
				new_mat[x][y] = matrix[x_name][y_name]

		return new_mat

if __name__ == "__main__":
	
	parser = ap.ArgumentParser()
	parser.add_argument('-d','--dist_mat',required=True)
	parser.add_argument('-p','--prefix',default="PHYLOGENTIC_FIGURE")
	parser.add_argument('-f','--format',default='.svg',nargs='*')
	parser.add_argument('-o','--out',default="PHYLOGENTIC_FIGURES")

	args = parser.parse_args()

	dist_mat = args.dist_mat
	prefix = args.prefix
	out_format = args.format
	outdir = args.out

	if not isdir(outdir):
		mkdir(outdir)

	Dij = []
	names = False

	try:
		FILE = open(dist_mat,"r")
		for line in FILE:
			line = line.replace("\n",'')
			if not names:
				names = line.split("\t")
			else:
				values = [float(i) for i in line.split("\t")]
				Dij.append(values)
	except:
		die(f"[E]  Unable to read from {dist_mat}")


	Tree(Dij,names,prefix,out_format,outdir)
