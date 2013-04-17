#!/usr/bin/env python

import sys,os,glob


if __name__ == "__main__":
	pair1,pair2,out = sys.argv[1:4]
	f1 = open(pair1)
	f2 = open(pair2)
	f0 = open(out,'w')
	last = None
	while last != f1.tell():
		last = f1.tell()
		for _ in range(4):
			f0.write(f1.readline())
		for _ in range(4):
			f0.write(f2.readline())
	f1.close()
	f2.close()
	f0.close()
	os.system('rm '+pair1+' '+pair2)