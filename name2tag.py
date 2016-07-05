import sys
import os
import re
import math


usage = """\n<tab delimited file to convert>\n\n"""
argv = sys.argv
junk = argv.pop(0)

if len(sys.argv) == 2:
	tag_file = argv.pop(0)
	tree = argv.pop(0)
else:
	sys.exit(usage)


TAG = {}
tag = open(tag_file, 'r')
tag = tag.readlines()
for i in tag:
	i = i.replace("\n",'')
	i = i.split("\t")
	if len(i) > 1:
		TAG[i[1]] = i[0]

tree = open(tree, 'r')
tree = tree.readlines()
for ln in tree:
	ln = ln.replace('\n','')
	ln_data = ln.split("\t")
	for i in ln_data:
		if TAG.has_key(i):
			ln = ln.replace(i,TAG[i])
	print ln


