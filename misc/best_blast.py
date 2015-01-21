import glob,os
from operator import itemgetter
from collections import defaultdict

def gen(fp,n=1):
	f = open(fp)
	G = defaultdict(float)
	for line in f:
		ls = line.strip().split('\t')
		if len(ls) == 2:
			G[' '.join(ls[0].split()[:n])] += float(ls[1])
	f.close()
	s = sum(G.values())
	X = [(k,v/s,s) for k,v in G.items()]
	return sorted(X,key=itemgetter(1),reverse=True)[0]

def assem_len(fp):
	f = open(fp[:fp.rfind('/')] + '/SPADES/scaffolds.fasta')
	l = 0
	for line in f:
		if line[:2] == '>N':
			l += int(line.split('_')[3])
	f.close()
	return l

FP = glob.glob(os.path.join('read_partitions/*/','alignments.bestHits.nameSum.txt'))
G = []
for fp in FP:
	g = gen(fp,n=2)
	n = fp[fp.index('/')+1:fp.index('/align')]
	G.append((g[0],g[1],g[2],assem_len(fp),n))

for g in sorted(G,key=itemgetter(2),reverse=True)[:10]:
	print g[0],g[1],g[2],g[3],g[4]
