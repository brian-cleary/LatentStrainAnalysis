import glob,os
from operator import itemgetter
from collections import defaultdict

def gen(fp):
	f = open(fp)
	G = defaultdict(float)
	for line in f:
		ls = line.strip().split('\t')
		if len(ls) == 2:
			G[ls[0].split()[0]] += float(ls[1])
	f.close()
	s = sum(G.values())
	X = [(k,v/s,s) for k,v in G.items()]
	return sorted(X,key=itemgetter(1),reverse=True)[0]

FP = glob.glob(os.path.join('read_partitions/*/','alignments.bestHits.nameSum.txt'))
G = []
for fp in FP:
	G.append(gen(fp))

for g in sorted(G,key=itemgetter(1),reverse=True)[:20]:
	print g[0],g[1],g[2]