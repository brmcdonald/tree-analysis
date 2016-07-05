import sys
import os
import re
from Bio import Phylo
sys.path.append('/home/bradon/scripts/BRM_LIB')
import brm_stdlib as brm
from scipy import stats

####################################################################################################
usage = """\n

Runs fisher's exact test and identifies gene families over/under represented in each subclade
of a phylogeny. Output is a matrix of gene families: 1 is uncorrelated, 2 is over represented,
3 is under represented.

<tree newick file>
<gene family table>
<max % of genomes to consider>
<min number of genomes for a subclade>
<max clusters to return, or ALL>
<output file>\n\n"""
argv = sys.argv[1:]

if len(argv) == 6:
	inTree = argv.pop(0)
	inFam = argv.pop(0)
	maxGenomes = float(argv.pop(0)) / 100
	minGenomes = float(argv.pop(0))
	cluster_limit = int(argv.pop(0))
	outFile = argv.pop(0)
else:
	sys.exit(usage)

class Clade(object):
	def __init__(self):
		self.a = {}
		self.b = {}
		self.c = {}
		self.d = {}
		self.mem = {}
		self.obj = ''

CLADES = {}
CLADE_GENE = {}
FAM_TOTAL = {}
TBL = brm.ortho2dict(inFam)

tree = Phylo.read(inTree,'newick')
totalTaxa = len(tree.get_terminals())

nonTerms = tree.get_nonterminals()
cladeCount = 0
for i in nonTerms:
	if len(i.get_terminals()) > minGenomes:
		x = Clade()
		CLADES[cladeCount] = x
		for g in i.get_terminals():
			g = g.name
			CLADES[cladeCount].mem[g] = {}
		CLADES[cladeCount].obj = i
		CLADES[cladeCount].obj.name = str(cladeCount)
	cladeCount += 1

for fam in TBL:
	if not FAM_TOTAL.has_key(fam):
		FAM_TOTAL[fam] = len(TBL[fam])
	for c in CLADES:
		CLADES[c].a[fam] = 0
		CLADES[c].b[fam] = 0
		CLADES[c].c[fam] = 0
		CLADES[c].d[fam] = 0

		if not CLADE_GENE.has_key(c):
			CLADE_GENE[c] = 0

		for g in TBL[fam]:
			genome = g.split('|')[0]

			#In clade, in fam
			if CLADES[c].mem.has_key(genome):
				CLADES[c].a[fam] += 1
				CLADE_GENE[c] += 1
			#Out clade, in fam
			else:
				CLADES[c].b[fam] += 1

for c in CLADES:
	for fam in CLADES[c].a:
		CLADES[c].c[fam] = FAM_TOTAL[fam] - CLADES[c].a[fam]
		CLADES[c].d[fam] = FAM_TOTAL[fam]  - CLADES[c].b[fam]

#Sanity check
for c in CLADES:
	total = 0
	for f in TBL:
		total += CLADES[c].a[f]
	if CLADE_GENE[c] != total:
		print c

OVER = {}
UNDER = {}
for clade in CLADES:
	OVER[clade] = {}
	UNDER[clade] = {}

#Running FET and picking out over/under represented clade/fam pairs
for clade in CLADES:
	perc_genomes = float(len(CLADES[clade].mem)) / float(totalTaxa)
	if perc_genomes <= maxGenomes and len(CLADES[clade].mem) > minGenomes:
		for fam in CLADES[c].a:
			a = CLADES[clade].a[fam]
			b = CLADES[clade].b[fam]
			c = CLADES[clade].c[fam]
			d = CLADES[clade].d[fam]
			oddsratio, pvalue = stats.fisher_exact([[a, b], [c, d]])

			if pvalue < 0.01 and oddsratio > 1:
				OVER[clade][fam] = 1
			elif pvalue < 0.01 and oddsratio < 1:
				UNDER[clade][fam] = 1

#Compiling matrix data from gene family table (modified from orth_tbl2matrix.py)
GENOME = {}
CAT = {}
ortho_data = open(inFam,'r')
tbl = ortho_data.readlines()
ortho_data.close()
for ln in tbl:
	ln = ln.split('\t')
	cat = ln[0].replace(':','')
	if not CAT.has_key(cat):
		CAT[cat] = {}

	genes = ln[1].replace('\n','').split(' ')
	for g in genes:
		if '|' in g:
			genome = g.split('|')[0]
		else:
			genome = g.split('_')[0]

		if not GENOME.has_key(genome):
			GENOME[genome] = {}

		if GENOME.has_key(genome):
			if not GENOME[genome].has_key(cat):
				GENOME[genome][cat] = 0
			GENOME[genome][cat] = 1
			CAT[cat][g] = ''

for g in GENOME:
	for cat in CAT:
		if not GENOME[g].has_key(cat):
			GENOME[g][cat] = 0

if cluster_limit == 'ALL':
	increment = 0
else:
	increment = 1

out = open(outFile,'w')
out.write('TBL')

for c in sorted(CLADES, key=lambda k: len(CLADES[k].mem.keys())):
	for cat in sorted(CAT, key=lambda k: len(CAT[k]), reverse=True):
		if OVER[c].has_key(cat):
			for g in CLADES[c].mem:
				if GENOME[g][cat] > 0:
					GENOME[g][cat] = 2
		elif UNDER[c].has_key(cat):
			for g in CLADES[c].mem:
				if GENOME[g][cat] > 0:
					pass
					#GENOME[g][cat] = 1

count = 0
for i in sorted(CAT, key=lambda k: len(CAT[k]), reverse=True):
	if count < cluster_limit:
		out.write('\t{0}'.format(i))
		count += increment
out.write('\n')


for g in GENOME:
	out.write(g)
	count = 0
	for c in sorted(CAT, key=lambda k: len(CAT[k]), reverse=True):
		if count < cluster_limit:
			out.write('\t{0}'.format(GENOME[g][c]))
			count += increment
	out.write('\n')
out.close()
