import sys,os,re
from multiprocessing import Process, Queue, Pool

usage = """
<folder of alignments>
(N)ucleotide or (P)rotein?>
<Num bootstraps (0 for no bootstrapping)>
<output folder>
<cores to use>
<cores per tree>\n\n"""
argv = sys.argv
junk = argv.pop(0)

if len(sys.argv) == 6:
	in_folder = argv.pop(0)
	type_flag = argv.pop(0)
	s_flag = int(argv.pop(0))
	out_folder = argv.pop(0)
	cores=int(argv.pop(0))
	coresPerTree = int(argv.pop(0))
else:
	sys.exit(usage)

def multi_run(cmd):
	os.system(cmd)

rxml = '~/tools_and_software/RAxML-8.1.24/raxmlHPC-PTHREADS-AVX2'

os.system('mkdir -p {0}'.format(out_folder))
aligns = os.listdir(in_folder)

cmd = []
NAMES = {}

for f in aligns:
	f_name = f
	f_name = f_name.replace('.fna','')
	f_name = f_name.replace('.faa','')
	f_name = f_name.replace('.mafft','')
	NAMES[f_name] = ''
	if type_flag == 'N':
		if s_flag > 0:
			data = (rxml,in_folder,f,f_name,coresPerTree,s_flag)
			cmd.append('{0} -s {1}/{2} -n {3} -x 3524627 -p 3524627 -f a -T {4} -N {5} -m GTRCAT -c 6'.format(*data))
		else:
			data = (rxml,in_folder,f,f_name,coresPerTree,s_flag)
			cmd.append('{0} -s {1}/{2} -n {3} -T {4} -m GTRCAT -c 6'.format(*data))
	elif type_flag == 'P':
		if s_flag > 0:
			data = (rxml,in_folder,f,f_name,coresPerTree,s_flag)
			cmd.append('{0} -s {1}/{2} -n {3} -x 3524627 -f a -T {4} -N {5} -m PROTGAMMABLOSUM62'.format(*data))
		else:
			data = (rxml,in_folder,f,f_name,coresPerTree,s_flag)
			cmd.append('{0} -s {1}/{2} -n {3} -T {4} -m PROTGAMMABLOSUM62'.format(*data))

numRun = cores / coresPerTree
p=Pool(numRun)
for c in cmd:
	p.apply_async(multi_run, args=(c,))
	pass
p.close()
p.join()

outputs = os.listdir('./')
for f in outputs:
	front = f.split('.')[0]
	tail = '.'.join(f.split('.')[1:])
	if s_flag == 0:
		if front == 'RAxML_bestTree' and NAMES.has_key(tail):
			os.system('mv {0} {1}/{2}.nwk'.format(f,out_folder,tail))
	else:
		if front == 'RAxML_bipartitionsBranchLabels' and NAMES.has_key(tail):
			os.system('mv {0} {1}/{2}.nwk'.format(f,out_folder,tail))

