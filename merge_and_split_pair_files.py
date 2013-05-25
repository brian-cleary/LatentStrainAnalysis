#!/usr/bin/env python

import sys,os,glob

def merge_pairs(f1,f2,f0):
	for i in range(0,25*10**6,10**5):
		r1 = [f1.readline() for _ in range(10**5)]
		r2 = [f2.readline() for _ in range(10**5)]
		if len(r1) > 0:
			while (r1[0].strip().split()[0][:-1] != r2[0].strip().split()[0][:-1]) and (r1[0][0] != '@'):
				r1 = r1[1:]
				r2 = r2[1:]
			for j in range(0,len(r1),4):
				f0.writelines(r1[j:j+4])
				f0.writelines(r2[j:j+4])
		else:
			break
	return len(r1)


if __name__ == "__main__":
	pair1,pair2,out = sys.argv[1:4]
	f1 = open(pair1)
	f2 = open(pair2)
	split_suffix = ['.0'+str(i) for i in range(10)]
	split_suffix += ['.'+str(i) for i in range(10,999)]
	r1len = 1
	ss = 0
	while r1len > 0:
		f0 = open(out+split_suffix[ss],'w')
		r1len = merge_pairs(f1,f2,f0)
		ss += 1
		f0.close()
	os.system('touch '+out)