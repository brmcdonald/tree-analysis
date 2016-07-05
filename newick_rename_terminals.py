import sys,os,re,math
from Bio import Phylo

usage = """

Changes labels of taxa in a tree. CANNOT deal with fasttree support scores. 
The final names cannot have single quotes in them.

<Lookup table>
<input tree>
<output tree>\n\n"""

try:
	argv = sys.argv[1:]
	inTable = argv.pop(0)
	inTree = argv.pop(0)
	outFile = argv.pop(0)
except IndexError:
	sys.exit(usage)
label = 0

tree = Phylo.read(inTree,'newick')

NAME = {}
with open(inTable,'r') as f:
	for ln in f:
		ln = ln.replace('\n','')
		ln = ln.split('\t')
		NAME[ln[0]] = ln[1]

terms = tree.get_terminals()
for i in terms:
	if NAME.has_key(i.name):
		i.name = NAME[i.name]


Phylo.write(tree,'zz_temp_rename_output.out','newick')
out = open(outFile,'w')
with open('zz_temp_rename_output.out','r') as f:
	for ln in f:
		ln = ln.replace('\'','')
		out.write(ln)
out.close()
os.system('rm zz_temp_rename_output.out')