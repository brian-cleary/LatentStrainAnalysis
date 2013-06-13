#!/usr/bin/env python

import sys,os,glob

def merge_pairs(f1,f2,f0):
	for i in range(0,8*10**6,10**5):
		r1 = [f1.readline() for _ in range(10**5)]
		r2 = [f2.readline() for _ in range(10**5)]
		try:
			while (r1[0].strip().split()[0][:-1] != r2[0].strip().split()[0][:-1]) and (r1[0][0] != '@'):
				r1 = r1[1:]
				r2 = r2[1:]
		except:
			r1 = ['']
			r2 = ['']
			pass
		for j in range(0,len(r1),4):
			f0.writelines(r1[j:j+4])
			f0.writelines(r2[j:j+4])
	return len(r1[0])

def split_singletons(sing_path,out_prefix):
	split_suffix = ['.0'+str(_) for _ in range(10)]
	split_suffix += ['.'+str(_) for _ in range(10,999)]
	ss = 0
	i = 0
	f1 = open(sing_path)
	for line in f1:
		if i%16000000 == 0:
			f0 = open(out_prefix+'.singleton.fastq'+split_suffix[ss],'w')
			ss += 1
		f0.write(line)
		i += 1


# python merge_and_split_pair_files.py A.fastq.1 A.fastq.2 A.fastq.sing original_reads/A
if __name__ == "__main__":
	pair1,pair2,sing,out = sys.argv[1:5]
	f1 = open(pair1)
	f2 = open(pair2)
	split_suffix = ['.0'+str(_) for _ in range(10)]
	split_suffix += ['.'+str(_) for _ in range(10,999)]
	r1len = 1
	ss = 0
	while r1len > 0:
		f0 = open(out+'.interleaved.fastq'+split_suffix[ss],'w')
		r1len = merge_pairs(f1,f2,f0)
		ss += 1
		f0.close()
	split_singletons(sing,out)
	os.system('touch '+out)