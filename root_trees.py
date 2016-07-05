import sys
import os
from Bio import Phylo


####################################################################################################
usage = """\n
Root a newick tree. Can either use a single terminal node as root or the common ancestor
of two terminal nodes, comma separated

<Tree file or folder of tree files>
<root location>
<output file/folder>\n\n"""
argv = sys.argv[1:]

if len(argv) == 3:
	inputData = argv.pop(0)
	rootTaxa = argv.pop(0)
	outLoc = argv.pop(0)
else:
	sys.exit(usage)

def rootTree(f, root,output):
	tree = Phylo.read(f,'newick')
	if ',' in root:
		taxa = root.split(',')
		root = tree.common_ancestor(taxa)
	tree.root_with_outgroup(root)
	Phylo.write(tree,output,'newick')

if os.path.isdir(inputData):
	trees = os.listdir(inputData)
	os.system('mkdir -p {0}'.format(outLoc))
	for t in trees:
		name = ''.join(t.split('.')[:-1])
		rootTree('{0}/{1}'.format(inputData,t),rootTaxa,'{0}/{1}_rooted.nwk'.format(outLoc,name))
else:
	name = inputData.split('/')[-1]
	print name
	rootTree('{0}'.format(inputData),rootTaxa,outLoc)

