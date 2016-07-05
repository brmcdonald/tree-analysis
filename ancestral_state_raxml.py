import sys
import os
import collections as cl
from Bio import Phylo

####################################################################################################
usage = """\n

Estimates ancestral states of internal nodes on a tree using RAxML. Cannot handle missing
taxa.

<Tree>
<Folder of alignments>
<List of alignments to process, or ALL>
<number of cores>
<output folder>\n\n"""

try:
	argv = sys.argv[1:]
	inTree = argv.pop(0)
	alignsFolder = argv.pop(0)
	alignList = argv.pop(0)
	coreNum = int(argv.pop(0))
	outFolder = argv.pop(0)
except IndexError:
	sys.exit(usage)


def fasta2dict(f, length_only=False):
	"""Generates a dictionary from a fasta file: dict[name] = sequence or dict[name]=seq_length"""

	DATA = {}
	if os.path.isdir(f):
		print "ERROR: cannot read a directory"
		return
	if length_only == False:
		with open(f,'r') as fasta:
			for ln in fasta:
				ln = ln.replace('\n','')
				if '>' in ln:
					name = ln.split()[0].replace('>','')
					DATA[name] = cl.deque()
				else:
					DATA[name].append(ln)
		for i in DATA:
			DATA[i] = ''.join(DATA[i])
	else:
		with open(f,'r') as fasta:
			for ln in fasta:
				ln = ln.replace('\n','')
				if '>' in ln:
					name = ln.split()[0].replace('>','')
					DATA[name] =0
				else:
					DATA[name] += len(ln)
		for i in DATA:
			DATA[i] = ''.join(DATA[i])
	return DATA


####################################################################################################

raxmlPath = "~/tools_and_software/RAxML-8.1.24/raxmlHPC-PTHREADS-SSE3"
os.system('mkdir -p {0}'.format(outFolder))
os.system('mkdir -p {0}/ancestral_states'.format(outFolder))
os.system('mkdir -p {0}/labeled_trees'.format(outFolder))

tree = Phylo.read(inTree,'newick')
fuckBiopython = tree.get_terminals()
terms = []
for i in fuckBiopython:
	terms.append(i.name)


GOOD = {}
if not alignList == 'ALL':
	with open(alignList,'r') as f:
		for ln in f:
			ln = ln.replace('\n','')
			if ln:
				ln = ln.split()[0]
				GOOD[ln] = ''

ALIGN = {}
aligns = os.listdir(alignsFolder)
for al in aligns:
	name = al.replace('W_','')
	name = name.replace('_align','')
	name = name.replace('.fasta','')
	name = name.replace('.mafft','')

	if alignList == "ALL":
		ALIGN[name] = "{0}/{1}".format(alignsFolder,al)
	elif GOOD.has_key(name):
		ALIGN[name] = "{0}/{1}".format(alignsFolder,al)

for i in ALIGN:
	SEQ = fasta2dict(ALIGN[i])
	length = len(SEQ[SEQ.keys()[0]])
	for ter in terms:
		if not SEQ.has_key(ter):
			SEQ[ter] = 'X' * length

	inputName = 'RAxML_inputSeq.{0}'.format(i)
	tempSeq = open(inputName,'w')
	for s in SEQ:
			tempSeq.write('>{0}\n{1}\n'.format(s,SEQ[s]))
	tempSeq.close()

	os.system('{0} -f A -t {1} -s {2} -m GTRGAMMA -T {4} -n {3}'.format(raxmlPath,inTree,inputName,i,coreNum))
	os.system('mv *marginalAncestralStates.{0} {1}/ancestral_states/{0}.ancestral_states'.format(i,outFolder))
	os.system('mv *nodeLabelledRootedTree.{0} {1}/labeled_trees/{0}.labeled_tree'.format(i,outFolder))
	os.system('rm RAxML*.{0}'.format(i))
	os.system('rm RAxML*.reduced')


