import glob,os
from operator import itemgetter
from numpy import argmax
from collections import defaultdict
from LSA import LSA

class Cluster_Analysis(LSA):

	def __init__(self,inputpath,outputpath,assembly_suffix='_velvet/',assembly_norm=2500000):
		super(Cluster_Analysis,self).__init__(inputpath,outputpath)
		self.assembly_suffix = assembly_suffix
		self.assembly_norm = assembly_norm
		# DUMB:
		self.read_length = 100

	def process_cluster(self,cluster_prefix):
		total_reads = 0
		singletons = []
		pair_channels = [[] for _ in range(self.num_libs)]
		for sample,lib in self.sample_library.items():
			try:
				sample_reads = self.sort_read_pairs(cluster_prefix+'.'+sample)
			except:
				sample_reads = 0
			if sample_reads > 0:
				singletons.append(self.input_path+'.'.join([cluster_prefix,sample,'singleton','fastq']))
				pair_channels[lib].append(self.input_path+'.'.join([cluster_prefix,sample,'pairs','fastq']))
				total_reads += sample_reads
		pair_channels = [pc for pc in pair_channels if len(pc)>0]
		ec = total_reads*self.read_length/self.assembly_norm
		self.assembly_loc = self.output_path+cluster_prefix+self.assembly_suffix
		self.assemble_cluster(cluster_prefix,singletons,pair_channels,exp_cov=max(ec,10),cov_cutoff=max(ec/15,2))
		h = []
		for i in range(29,61,10):
			try:
				L = open(self.assembly_loc+'_'+str(i)+'/Log').readlines()
				h.append(int(L[-1].split(',')[0].split()[-1]))
			except:
				h.append(0)
		amax = argmax(h)
		os.system('cp '+self.assembly_loc+'_'+str(range(29,61,10)[amax])+'/contigs.fa '+self.assembly_loc+'contigs.fa')
		self.align_assembly(cluster_prefix)
		self.phyler_assembly(cluster_prefix)

	def assemble_cluster(self,cluster_prefix,sg,pc,exp_cov=50,cov_cutoff=15,min_contig=500):
		os.system('mkdir '+self.assembly_loc)
		i = 1
		s = []
		for c in pc:
			s.append('-shortPaired')
			if i > 1:
				s[-1] += str(i)
			s += c
			i += 1
		os.system('/home/unix/bcleary/src/velvet/velveth '+self.assembly_loc+' 29,61,10 -fastq -short '+' '.join(sg)+' '+' '.join(s))
		for i in range(29,61,10):
			os.system('/home/unix/bcleary/src/velvet/velvetg '+self.assembly_loc+'_'+str(i)+'/ -exp_cov '+str(exp_cov)+' -cov_cutoff '+str(cov_cutoff)+' -min_contig_lgth '+str(min_contig))

	def assemble_cluster_idba(self,cluster_prefix,min_contig=500):
		read_loc = self.input_path+cluster_prefix
		os.system('/bli/bin/Linux/x86_64/idba_ud-1.0.9/bin/fq2fa '+read_loc+'.pairs.fastq '+read_loc+'.pairs.fasta')
		os.system('mkdir '+self.assembly_loc)
		os.system('/bli/bin/Linux/x86_64/idba_ud-1.0.9/bin/idba_ud -r '+read_loc+'.pairs.fasta -o '+self.assembly_loc+' --min_contig '+str(min_contig)+' --pre_correction')

	def align_assembly(self,cluster_prefix):
		os.system('/bli/bin/Linux/x86_64/ncbi-blast-2.2.27/ncbi-blast-2.2.27+/bin/blastn -query '+self.assembly_loc+'contigs.fa -db /import/pool2/projects/mega/search_dbs/prod/BLIS_SEQUENCES/blast/Bacterial_DNA_All_Prokaryotic -task megablast -outfmt "7 qseqid qstart qend sseqid score" -out '+self.output_path+cluster_prefix+'.alignments.txt')
		os.system('/bli/bin/Linux/x86_64/ncbi-blast-2.2.27/ncbi-blast-2.2.27+/bin/blastn -query '+self.assembly_loc+'contigs.fa -db /import/pool2/projects/mega/search_dbs/prod/BLIS_SEQUENCES/blast/Viral_DNA_All_GenBank_Virus_Sequences -task megablast -outfmt "7 qseqid qstart qend sseqid score" -out '+self.output_path+cluster_prefix+'.viral.alignments.txt')

	def phyler_assembly(self,cluster_prefix):
		os.system('/home/unix/bcleary/src/MetaPhylerV1.25/runMetaphyler.pl '+self.assembly_loc+'contigs.fa blastn '+self.output_path+cluster_prefix+'.phyler.blastn 4')

	def assembly_stats(self,cluster_prefix,file_suffix='contigs.fa'):
		AS = {}
		L = []
		f = open(self.assembly_loc+file_suffix)
		line = f.readline()
		while line != '':
			if line[0] == '>':
				L.append(int(line.replace('_',' ').split()[3]))
			line = f.readline()
		f.close()
		if L:
			L.sort(reverse=True)
			t = 0
			for l in L:
				t += l
				if t > sum(L)/2:
					break
			AS = {'total bp': sum(L),'num scaffolds': len(L),'N50': l,'largest contig': L[0]}
		return AS

	def alignment_stats(self,cluster_prefix):
		f = open(self.output_path+cluster_prefix+'.alignments.txt')
		L = f.readlines()
		f.close()
		A = defaultdict(int)
		current_alignments = {}
		for l in L:
			lst = l.strip().split('\t')
			if len(lst) != 5:
				for k,v in current_alignments.items():
					A[k[0]] += v
				current_alignments = {}
			else:
				alignid_querypos = (lst[3],(lst[1],lst[2]))
				if alignid_querypos not in current_alignments:
					current_alignments[alignid_querypos] = int(lst[4])
		return A

	def phyler_stats(self,cluster_prefix,level='genus',type='blastn'):
		L = open(self.output_path+cluster_prefix+'.'.join(['.phyler',type,level,'taxprof'])).readlines()
		A = []
		for l in L[1:]:
			x = l.strip().split()
			A.append((x[0],float(x[1]),int(x[2])))
		return A

	def sort_read_pairs(self,cluster_prefix):
		try:
			f = open(self.input_path+cluster_prefix+'.fastq','r')
		except:
			return 0
		type = self.id_type(f)
		# this is kind of a bummer since files are *mostly* sorted already
		if type == 1:
			sorted_reads = sorted(self.read_generator(f,raw_reads=True,max_reads=10**7),key=lambda (d): d[:d.index(' ')+2])
			def get_id(r):
				return r[:r.index(' ')+2]
		elif type == 2:
			sorted_reads = sorted(self.read_generator(f,raw_reads=True,max_reads=10**7),key=lambda (d): d.split()[0])
			def get_id(r):
				return r.split()[0]
		else:
			sorted_reads = sorted(self.read_generator(f,raw_reads=True,max_reads=10**7),key=lambda (d): d.split()[0])
			def get_id(r):
				return r.split()[0]+'*'
		if len(sorted_reads) > 0:
			pair_file1 = open(self.input_path+cluster_prefix+'.mate1.fastq','w')
			pair_file2 = open(self.input_path+cluster_prefix+'.mate2.fastq','w')
			singleton_file = open(self.input_path+cluster_prefix+'.singleton.fastq','w')
			total_reads = 0
			while len(sorted_reads) > 0:
				current_id = ''
				last = ''
				for r in sorted_reads:
					r_id = get_id(r)
					if (r_id[:-1] == current_id[:-1]) and (r_id != current_id):
						pair_file1.write(last)
						pair_file2.write(r)
						total_reads += 2
						last = ''
					else:
						singleton_file.write(last)
						total_reads += 1
						last = r
					current_id = r_id
				if type == 1:
					sorted_reads = sorted(self.read_generator(f,raw_reads=True,max_reads=10**7),key=lambda (d): d[:d.index(' ')+2])
				elif type == 2:
					sorted_reads = sorted(self.read_generator(f,raw_reads=True,max_reads=10**7),key=lambda (d): d.split()[0])
				else:
					sorted_reads = sorted(self.read_generator(f,raw_reads=True,max_reads=10**7),key=lambda (d): d.split()[0])
			pair_file1.close()
			pair_file2.close()
			singleton_file.close()
			f.close()
			os.system('rm '+self.input_path+cluster_prefix+'.fastq')
			os.system('touch '+self.input_path+cluster_prefix+'.fastq')
			return total_reads
		else:
			f.close()
			try:
				r = sum(self.read_counts(cluster_prefix).values())
			except:
				r = 0
			return r

	def read_counts(self,cluster_prefix):
		RC = {}
		RC['paired'] = self.read_count(self.input_path+cluster_prefix+'.mate1.fastq',startchar='@')*2
		RC['singleton'] = self.read_count(self.input_path+cluster_prefix+'.singleton.fastq',startchar='@')
		return RC

	def read_count(self,fp,startchar='>'):
		f = open(fp,'r')
		r = 0
		line = f.readline()
		while line != '':
			if line[0] == startchar:
				r += 1
			line = f.readline()
		f.close()
		return r