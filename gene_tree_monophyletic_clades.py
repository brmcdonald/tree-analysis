import sys
import os
from Bio import Phylo


####################################################################################################
usage = """\n
Given a set of trees, identify trees in which one of the given groups of strains is monophyletic.
Can deal with unrooted trees with missing taxa. Only works on unrooted trees with at least 3 
internal nodes, otherwise midpoint rooting fails.

Clade file should be fasta-esqe, ie:

>Clade1
strain1
strain2
>Clade2
strain3
strain4

<Folder of gene trees>
<clade file>
<Unrooted trees? (Y/N)>
<Minimum support value>
<Annotation table for gene names>
<Min organisms for a monophyletic clade>
<Min genes in a cluster>
<Max distance between genes in a cluster>
<output folder>\n\n"""

try:
	argv = sys.argv[1:]
	inputData = argv.pop(0)
	inClade = argv.pop(0)
	rootFlag = argv.pop(0).upper()
	minConf = float(argv.pop(0))
	annTable = argv.pop(0)
	minGenomes = int(argv.pop(0))
	minCluster = int(argv.pop(0))
	maxDist = int(argv.pop(0))
	outLoc = argv.pop(0)
except IndexError:
	sys.exit(usage)

class Clade(object):
	def __init__(self):
		self.members = {}		#Which genomes are found in the clade
		self.monoTrees = {}	#Which trees meet monophyly criteria

		self.terminals = [] #Get terminal node objects because biopython is shit

	def resetTerm(self):
		self.terminals = []

class Cluster(object):
	def __init__(self):
		self.genes = []
		self.last = 0

os.system('mkdir -p {0}'.format(outLoc))

#Reading annotation table
ANN = {}
with open(annTable,'r') as f:
	for ln in f:
		ln = ln.split('\t')
		name = ln[0].split('|')[-1]
		gene = ln[1]
		ANN[name] = gene

#Reading clades
CLADE = {}
with open(inClade,'r') as f:
	for ln in f:
		ln = ln.replace('\n','')
		if '>' in ln:
			clade = ln.replace('>','').split()[0]
			CLADE[clade] = Clade()
		elif ln:
			CLADE[clade].members[ln] = ''

trees = os.listdir(inputData)
for t in trees:
	n = t.replace('.fasttree','')
	n = t.replace('.fastttree','')
	n = n.replace('_align','')
	n = n.replace('W_','')
	n = n.replace('.nwk','')
	n = n.replace('.newick','')

	nwk = Phylo.read('{0}/{1}'.format(inputData,t),'newick')
	if len(nwk.get_nonterminals()) >= 3:
		if rootFlag == 'Y':
			nwk.root_at_midpoint()

		for c in CLADE:
			#Workaround for Biopython being shit and forceing me to
			#search for node objects rather than just names
			CLADE[c].resetTerm()
			terms = nwk.get_terminals()
			for i in terms:
				if CLADE[c].members.has_key(i.name):
					CLADE[c].terminals.append(i)

			#Now I can actually search nodes
			try:
				if nwk.is_monophyletic(CLADE[c].terminals) and len(CLADE[c].terminals) >= minGenomes:
					node = nwk.common_ancestor(CLADE[c].terminals)
					if not node.confidence:
						CLADE[c].monoTrees[n] = 'Poly'
					elif node.confidence >= minConf:
						CLADE[c].monoTrees[n] = node.confidence
			except ValueError:
				pass

#Printing genes with monophyletic trees
for c in CLADE:
	out = open('{0}/{1}.monoTrees'.format(outLoc,c),'w')
	for t in sorted(CLADE[c].monoTrees):
		out.write('{0}\t{1}\t{2}\n'.format(t,ANN[t],CLADE[c].monoTrees[t]))
	out.close()

#Finding clusters of genes
raw = os.listdir(outLoc)
os.system('mkdir -p {0}/clusters'.format(outLoc))
for dat in raw:
	with open('{0}/{1}'.format(outLoc,dat)) as f:
		name = dat.replace('.monoTrees','')
		DATA = []
		last = 0
		newC = []
		for ln in f:
			ln = ln.replace('\n','')
			gData = ln
			ln = ln.split()
			gName = ln[0]
			gFunc = ln[1]
			gNum = int(gName.split('_')[-1])

			if last > 0:
				if gNum - last <= maxDist:
						newC.append(gData)
				else:
					if len(newC) >= minCluster:
						DATA.append(newC)
					newC = [gData]
			last = gNum
	out = open('{0}/clusters/{1}.clusters'.format(outLoc,name),'w')
	for c in DATA:
		out.write('>\n')
		for i in c:
			out.write(i + '\n')
	out.close()