import sys,os,re,math

usage = """\nLabels the internal nodes of a newick tree with temporary labels
<input tree>
<output tree>\n\n"""
argv = sys.argv
junk = argv.pop(0)

if len(sys.argv) == 2:
	inputTree = argv.pop(0)
	outTree = argv.pop(0)
else:
	sys.exit(usage)

label = 0

inTree = open(inputTree, 'r')
inTree = inTree.readlines()
output = open(outTree, "w")


for ln in inTree:
	ln = ln.split("):")
	tail = ln.pop()
	for i in ln:
		output.write("{0})A{1}:".format(i,label))
		label += 1
	output.write(tail)
output.close()
