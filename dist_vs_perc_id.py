import sys
import os
import ete2
#from Bio import Phylo

####################################################################################################
usage = """\n


<Percent identity table>
<phylogenetic tree>
<root taxa of tree>
<(T)ree distance or (D)ivergence time?> 
<output file>\n\n"""
argv = sys.argv[1:]

if len(argv) == 5:
	inFile = argv.pop(0)
	treeFile = argv.pop(0)
	rootTaxa = argv.pop(0)
	mode = argv.pop(0).upper()
	outFile = argv.pop(0)
else:
	sys.exit(usage)


def pairDist(tree,taxa):
	DIST = {}

	for i in taxa:
		DIST[i] = []
		for j in taxa:
			if i != j:
				pairdist = tree.get_distance(i,j)
				DIST[i].append(pairdist)
	return DIST


tree = ete2.Tree(treeFile)
#tree.set_outgroup(rootTaxa)


ID = {}
with open(inFile,'r') as f:
	for ln in f:
		ln = ln.replace('\n','')
		ln = ln.split()
		if not ID.has_key(ln[0]):
			ID[ln[0]] = {}
		ID[ln[0]][ln[1]] = float(ln[2])

DIST = {}
for x in ID:
	DIST[x] = {}
	for y in ID[x]:
		if x != y:
			if mode == 'T':
				DIST[x][y] = tree.get_distance(x,y)
			elif mode == 'D':
				DIST[x][y] = (tree.get_distance(x,y) / 2)
out = open(outFile,'w')
for x in DIST:
	for y in DIST[x]:
		if x != y:
			out.write('{0}\t{1}\t{2}\t{3}\n'.format(x,y,ID[x][y],DIST[x][y]))
out.close()
