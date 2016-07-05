import sys
import os
import collections as cl
import ete2
from scipy import stats
import numpy as np
sys.path.append('/home/bradon/scripts/BRM_LIB')
import brm_stdlib as brm

####################################################################################################
usage = """\n
Permutation branch length test for non-random association between genes and tree taxa. Note
that it only tests presence/absence of a gene family in a genome.

The analysis can be limited to subsets of the tree by specificing which taxa should be
included. In this case only genes found in those taxa will be counted, and all permuted trees
will be generated from members of that set of taxa.

<newick tree>
<orthoMCL style tbl>
<Minimum number of genes to run test on>
<Maximum number of genes to run test on>
<Number of random tree permutations>
<name for output files>
<Optional file of taxa to include>\n\n"""
argv = sys.argv[1:]

if len(argv) == 6:
	treeFile = argv.pop(0)
	geneTbl = argv.pop(0)
	minGenes = int(argv.pop(0))
	maxGenes = int(argv.pop(0))
	permNum = int(argv.pop(0))
	outFile = argv.pop(0)
	taxaList = ''
elif len(argv) == 7:
	treeFile = argv.pop(0)
	geneTbl = argv.pop(0)
	minGenes = int(argv.pop(0))
	maxGenes = int(argv.pop(0))
	permNum = int(argv.pop(0))
	outFile = argv.pop(0)
	taxaList = argv.pop(0)
else:
	sys.exit(usage)

#Generates permutations of the tree with a set number of taxa and returns a list of total
#branch lengths for the permuted trees
def permTree(tree,num_taxa,num_perm,debug=0):

	global TAXA
	lengths = cl.deque()
	leaves = []

	if not TAXA.keys():
		crap_leaves = tree.get_leaves()
		for i in crap_leaves:
			leaves.append(i.name)
	elif TAXA.keys():
		leaves = TAXA.keys()

	for x in range(0,num_perm):
		perm_tree = tree.copy()
		np.random.shuffle(leaves)
		good_leaves = leaves[:num_taxa]
		perm_tree.prune(good_leaves,preserve_branch_length=True)

		total = 0
		all_branch = perm_tree.traverse()
		for i in all_branch:
			total += i.dist
		lengths.append(total)
	return lengths


if taxaList:
	TAXA = {}
	with open(taxaList,'r') as f:
		for ln in f:
			if '>' not in ln:
				ln = ln.replace('\n','')
				ln = ln.split()[0]
				TAXA[ln] = ''

#Counting the number of taxa for each family, only running the tree permutation for each number
#of taxa once
tCore = ete2.Tree(treeFile)
with open(geneTbl,'r') as f:
	needPerm = []
	fams = 0
	for ln in f:
		name = ln.split('\t')[0].replace(':','')
		genes = ln.split('\t')[1].replace('\n','').split(' ')
		genomes = set([x.split('|')[0] for x in genes])

		#removing unwanted taxa if required
		if taxaList:
			removeList = []
			for i in genomes:
				if not TAXA.has_key(i):
					removeList.append(i)
		for i in removeList:
			genomes.remove(i)

		if len(genomes) >= minGenes and len(genomes) <= maxGenes:
			fams += 1
			needPerm.append(len(genomes))
	needPerm = set(needPerm)

print '\nExecuting {2} permutations of {0} trees for {1} gene fams...'.format(len(needPerm),fams,permNum)
PERM = {}
progress = 0
for i in needPerm:
	PERM[i] = permTree(tCore,i,permNum)
	progress += 1
	if progress % 10 == 0:
		print '\t{0} trees done'.format(progress)

print 'Calculating T-Test...'
out = open(outFile + '_ttest.out','w')
out.write('Fam\tPval\tT_stat\tfam_len\tavg_rand_len\tRand_stdev\tNum_StDev\tNum_taxa\n')
with open(geneTbl,'r') as f:
	for ln in f:
		name = ln.split('\t')[0].replace(':','')
		genes = ln.split('\t')[1].replace('\n','').split(' ')
		genomes = set([x.split('|')[0] for x in genes])

		#removing unwanted taxa if required
		if taxaList:
			removeList = []
			for i in genomes:
				if not TAXA.has_key(i):
					removeList.append(i)
		for i in removeList:
			genomes.remove(i)

		if len(genomes) >= minGenes and len(genomes) <= maxGenes:

			#Getting random tree length distribution
			ranDist = PERM[len(genomes)]
			avgDist = np.mean(ranDist)
			stdevDist = np.std(ranDist)

			#Calculating actual distribution based on genome distances
			perm_tree = tCore.copy()
			perm_tree.prune(genomes,preserve_branch_length=True)
			total = 0
			all_branch = perm_tree.traverse()
			for i in all_branch:
				total += i.dist
			num_stdev = (total - avgDist) / stdevDist

			t, pval = stats.ttest_1samp(ranDist, total)
	
			printData = (name,pval,t,total,avgDist,stdevDist,num_stdev,len(genomes))
			out.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(*printData))
out.close()


