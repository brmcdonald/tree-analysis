import sys,os,re
from multiprocessing import Process, Queue, Pool

usage = """
<folder of alignments>
<(N)ucleotide or (P)rotein?>
<compute support values? (Y/N)
<output folder>
<cores to use>\n\n"""
argv = sys.argv
junk = argv.pop(0)

if len(sys.argv) == 5:
	in_folder = argv.pop(0)
	type_flag = argv.pop(0)
	s_flag = argv.pop(0)
	out_folder = argv.pop(0)
	cores=int(argv.pop(0))
else:
	sys.exit(usage)

def multi_run(cmd):
	os.system(cmd)

ft_nt = '~/tools_and_software/fasttree/FastTree -nt -gtr'
ft_p = '~/tools_and_software/fasttree/FastTree'

if s_flag == 'N':
	ft_nt += ' -nosupport'
	ft_p += ' -nosupport'

os.system('mkdir -p {0}'.format(out_folder))
aligns = os.listdir(in_folder)

cmd = []

for f in aligns:
	f_name = f
	f_name = f_name.replace('.fna','')
	f_name = f_name.replace('.faa','')
	f_name = f_name.replace('.mafft','')
	f_name = f_name.replace('.fasta','')

	if type_flag == 'N':
		cmd.append('{0} < ./{1}/{2} > {3}/{4}.fastttree'.format(ft_nt,in_folder,f,out_folder,f_name))
	elif type_flag == 'P':
		cmd.append('{0} < ./{1}/{2} > {3}/{4}.fasttree'.format(ft_p,in_folder,f,out_folder,f_name))

p=Pool(cores)
for c in cmd:
	p.apply_async(multi_run, args=(c,))
	pass
p.close()
p.join()
