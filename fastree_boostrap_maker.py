import sys
import os
import re
import math
from multiprocessing import Process, Queue, Pool

usage = """\nUses RAxML to generate resampled alignments, 
then uses FASTTREE to generate trees of each
\n<folder of alignments>
<sequence type: dna or aa>
<output folder>
<Num bootstraps>
<num cores>\n\n"""
argv = sys.argv
junk = argv.pop(0)

if len(sys.argv) == 5:
	folder = argv.pop(0)
	seqType = argv.pop(0).lower()
	output_folder = argv.pop(0)
	boot_num = int(argv.pop(0))
	cores = int(argv.pop(0))
else:
	sys.exit(usage)

#Environmental variables

#Fasta to phylip format script
f2p = "perl ~/scripts/fasta_manipulation/fasta_alignment2phylip.pl"
#Raxml path
rax = "~/tools_and_software/RAxML-8.1.24/raxmlHPC-PTHREADS-AVX2"

def build_tree(f,seqType):
	os.system("perl ~/scripts/fasta_manipulation/phylip2fasta.pl {0}/raw_data/{1} > {0}/raw_data/{1}.fasta".format(output_folder,f))
	os.system('rm -f {0}/raw_data/{1}'.format(output_folder,f))
	if seqType == 'dna':
		os.system("~/tools_and_software/fasttree/FastTree -nt -gtr -nosupport < ./{0}/raw_data/{1}.fasta > {0}/rough_trees/{1}".format(output_folder,f))
	elif seqType == 'aa':		
		os.system("~/tools_and_software/fasttree/FastTree -nosupport < ./{0}/raw_data/{1}.fasta > {0}/rough_trees/{1}".format(output_folder,f))

os.system("mkdir -p {0}/raw_data".format(output_folder))
os.system("mkdir -p {0}/rough_trees".format(output_folder))
os.system("mkdir -p {0}/boot_trees".format(output_folder))

files = os.listdir(folder)
total = len(files)
per_folder = (total / cores)

name_lst = []
for f in files:
	name = f
	name = name.replace(".fna","")
	name = name.replace(".fasta","")
	name = name.replace(".mafft","")
	name = name.replace(".faa","")

	if os.path.getsize('{0}/{1}'.format(folder,f)): 
		os.system("{0} {1}/{2} > {3}/raw_data/{4}.BS_0".format(f2p,folder,f,output_folder,name))

		call_data = (rax,cores,(boot_num-1),name,output_folder)
		if seqType == 'dna':
			os.system("{0} -f j -T {1} -b 9 -# {2} -m GTRGAMMA -n {3} -s {4}/raw_data/{3}.BS_0".format(*call_data))
		elif seqType == 'aa':
			os.system("{0} -f j -T {1} -b 9 -# {2} -m PROTGAMMABLOSUM62 -n {3} -s {4}/raw_data/{3}.BS_0".format(*call_data))
		#os.system("mv *{0}.B* {1}/raw_data/".format(name,output_folder))

		#Remove useless stuff
		os.system("rm -f RAxML_info*")
		os.system("rm -f {3}/raw_data/{4}.phylip".format(f2p,folder,f,output_folder,name))
		os.system("rm -f {3}/raw_data/*.reduced".format(f2p,folder,f,output_folder,name))

files = os.listdir("./{0}/raw_data".format(output_folder))
p=Pool(cores)
for i in files:
	if ".BS" in str(i) and 'reduced' not in str(i):
		name = i.split('.')
		name = name[0]
		name_lst.append(name)
		p.apply_async(build_tree, args=(i,seqType))
		pass
	pass
p.close()
p.join()

name_lst = set(name_lst)
all_boots = os.listdir('{0}/rough_trees/'.format(output_folder))
for n in name_lst:
	boots_1 = [f for f in all_boots if '{0}.'.format(n) in f]
	boots = ' '.join(['./{0}/rough_trees/{1}'.format(output_folder,f) for f in boots_1])
	os.system("cat {3} > {0}/boot_trees/{1}_boot_{2}".format(output_folder,n,boot_num,boots))
#os.system("rm -f -r {0}/raw_data".format(output_folder))
#os.system("rm -f -r {0}/rough_trees".format(output_folder))
