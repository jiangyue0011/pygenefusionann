#!/usr/bin/env python

import sys
import pysam
import bisect
import time
import sequtils
print >> sys.stderr, "Version 0.2"
class GeneBed():
	def __init__(self, bed_line):
		tmp = bed_line.split()
		self.chr = tmp[0]
		self.start = int(tmp[1]) + 1	# bed file use 0-based coordinate
		self.end = int(tmp[2])		# start and end are first and last base of each segment
		self.transcript_id = tmp[3]
		self.type = tmp[4]
		self.idx = int(tmp[5])
		self.strand = tmp[6]
		self.gene_name = tmp[7]
		self.gene_id = tmp[8]
	def __tostring__(self):
		return bed_line.strip()
class CffFileStats():
	__fusion_dict = {}
	__fusion_samples_dict = {}
	def __init__(self, cff_file):
		self.__load_fusions(cff_file)
	def __load_fusions(self, cff_file):
		for line in open(cff_file, "r"):
			fusion = CffFusion(line)
	
			key = fusion.gene1 + "_" + fusion.gene2
			rkey = fusion.gene2 + "_" + fusion.gene1
			if rkey not in self.__fusion_dict:
				self.__fusion_dict.setdefault(key, []).append(fusion)
				self.__fusion_samples_dict.setdefault(key, []).append(fusion.sample)
			else:
				self.__fusion_dict.setdefault(rkey,[]).append(fusion)
				self.__fusion_samples_dict.setdefault(rkey, []).append(fusion.sample)
	def get_gene_order_stats(self):
		for key in self.__fusion_dict:
			n_sg = 0
			n_rt = 0 # read through
			n_gf = 0 # gene fusion
			n_tc = 0 # truncated coding
			n_tn = 0 # truncated noncoding
			n_ns = 0 # nonsense
			gene_order = []
			sample_type = []
			for sample in set(self.__fusion_samples_dict[key]):
				if sample.endswith("N"):
					sample_type.append("Normal")
				elif sample.endswith("T"):
					sample_type.append("Tumor")
				else:
					sample_type.append("Unknown")

			type = "Unknown"
			for fusion in self.__fusion_dict[key]:
				tmp = ";".join(fusion.otherann)
				gene_order.append(tmp)
				if "SameGene" in tmp:
					type = "SameGene"
					n_sg += 1
				elif "ReadThrough" in tmp:
					type = "ReadThrough" if type == "Unknown" else type
					n_rt += 1
				elif "GeneFusion" in tmp:
					type = "GeneFusion" if type == "Unknown" else type
					n_gf += 1
				elif "TruncatedCoding" in tmp:
					type = "TruncatedCoding" if type == "Unknown" else type
					n_tc +=1
				elif "TruncatedNoncoding" in tmp:
					type = "TruncatedNoncoding" if type == "Unknown" else type
					n_tn += 1
				elif "NoDriverGene" in tmp:
					type = "NoDriverGene" if type == "Unknown" else type
					n_ns += 1
				else:
					print >> sys.stderr, "Fusions without category:", fusion.tostring()
			if type != "Unknown":
				print "Fusion", fusion.gene1, fusion.gene2, ",".join(list(set(self.__fusion_samples_dict[key]))), type, ",".join(sorted(list(set(sample_type)))),
				print "|".join(list(set(gene_order)))
				print "\tSameGene:", n_sg
				print "\tReadThrough:", n_rt
				print "\tGeneFusion:", n_gf
				print "\tTruncatedCoding:", n_tc
				print "\tTruncatedNoncoding:", n_tn
				print "\tNonSense:", n_ns


class CffFusion():
	def __init__(self, cff_line):
		tmp = cff_line.split()
		self.chr1 = tmp[0]
		self.pos1 = int(tmp[1])
		self.strand1 = tmp[2]
		self.chr2 = tmp[3]
		self.pos2 = int(tmp[4])
		self.strand2 = tmp[5]
		self.orf = tmp[6]
		self.read_through = tmp[7]
		self.split_cnt = int(tmp[8])
		self.span_cnt = int(tmp[9]) if tmp[9] != "NA" else tmp[9]
		self.sample = tmp[10]
		self.lib = tmp[11]
		self.tool = tmp[12]
		self.id = tmp[13]
		self.score = float(tmp[14]) if tmp[14] != "NA" else tmp[9]
		self.otherann = tmp[19:]
		self.gene1 = tmp[15]
		self.type1 = tmp[16]
		self.gene2 = tmp[17]
		self.type2 = tmp[18]
		self.trans_id1 = "NA"
		self.trans_id2 = "NA"
		#self.elements = tmp[0:]

		if not self.chr1.startswith("chr"):
			self.chr1 = "chr" + self.chr1
			self.chr2 = "chr" + self.chr2
	def tostring(self):
		#return " \t".join(self.elements) + "\t" +  self.gene1 + "\t" + self.gene2
		return "\t".join(map(lambda x:str(x), [self.chr1, self.pos1, self.strand1, self.chr2, self.pos2, self.strand2, self.orf, self.read_through, self.split_cnt, self.span_cnt, self.sample, self.lib, self.tool, self.id, self.score, self.gene1, self.type1, self.gene2, self.type2, "\t".join(self.otherann)]))
	# according to fusion strand (defuse style, strands are supporting pairs') return all possible gene fusions
	def __check_gene_pairs(self, genes1, genes2, gene_ann):
		gene_order = []
		type1 = []
		type2 = []
		id1 = []
		id2 = []
		category = ""
		common_genes =  set(genes1.keys()) & set(genes2.keys())
		for gene_name in set(genes1.keys()):
			if gene_ann.is_coding(gene_name):
				type1.append("CodingGene")
				id1.append(str(gene_ann.get_coding_gene_idx(gene_name)))
			else:
				type1.append("NoncodingGene")
		for gene_name in set(genes2.keys()):
			if gene_ann.is_coding(gene_name):
				type2.append("CodingGene")
				id2.append(str(gene_ann.get_coding_gene_idx(gene_name)))
			else:
				type2.append("NoncodingGene")

		# No driver gene
		if not genes1:
			gene_order.append("NoDriverGene")
			gene_order.extend(["NA"]*3)
		# map to same gene
		elif common_genes:
			gene_order.append("SameGene" + "\t" + ",".join(list(common_genes)))
			gene_order.extend(["NA"]*2)
		else:
			list1 = []
			for gene_name in genes1:
				tmp = set(genes1[gene_name])
				list1.append(gene_name + "_" + ",".join(list(tmp)))
			list2 = []
			for gene_name in genes2:
				tmp = set(genes2[gene_name])
				list2.append(gene_name + "_" + ",".join(list(tmp)))
			#gene_order.append(",".join(list(set(genes1))) + ">>" + ",".join(list(set(genes2))))

			# category fusions into: read through, gene fusion, truncated coding, truncated noncoding, nonsense
			if "CodingGene" in type1 and "CodingGene" in type2:
				for id in id1:
					tmp = id.split("_")
					id = int(tmp[0])
					strand = tmp[1]
					if (strand == "f" and  str(id+1) + "_f" in id2 ) or (strand == "r" and  str(id-1) + "_r" in id2):
						category = "ReadThrough"
					else:
						category = "GeneFusion"
			elif "CodingGene" in type1:
				category = "TruncatedCoding"
			elif "NoncodingGene" in type1:
				category = "TruncatedNoncoding"
			elif not type1:
				category = "NonSense"
			else:
				print >> sys.stderr, "Warning: Unknown category."
				print >> sys.stderr, type1, type2



			gene_order.append(",".join(list1) + ">>" + ",".join(list2))
			gene_order.append(",".join(type1) + ">>" + ",".join(type2))
			gene_order.append(",".join(id1) + ">>" + ",".join(id2))
			gene_order.append(category)

		return gene_order

	def ann_gene_order(self, gene_ann):
		gene_order = []
		# fusion has been annotated
		if self.otherann:
			return gene_order
		matched_genes1 = gene_ann.map_pos_to_genes(self.chr1, self.pos1)
		matched_genes2 = gene_ann.map_pos_to_genes(self.chr2, self.pos2)

		a = {} # forward strand gene at pos1
		c = {} # backward strand gene at pos1
		b = {} # forward strand gene at pos2
		d = {} # backward strand gene at pos2
		for gene in matched_genes1:
			if gene.strand == "f":
				a.setdefault(gene.gene_name, []).append(gene.type)
			else:
				c.setdefault(gene.gene_name, []).append(gene.type)
		for gene in matched_genes2:
			if gene.strand == "f":
				b.setdefault(gene.gene_name, []).append(gene.type)
			else:
				d.setdefault(gene.gene_name, []).append(gene.type)
		# for tools do not provide defuse-style strand, regenerate strands, this is assuming that gene1 = 5 prime gene and gene2 = 3 prime gene
		if self.strand1 == "NA" or self.strand2 == "NA":
			gene_interval1 = ""
			gene_interval2 = ""
			for sep in [",", "/"]:
				for gene_name in (self.gene1).split(sep):
					if not gene_interval1:
						gene_interval1 = gene_ann.get_gene_interval(gene_name)
					else:
						break
				for gene_name in (self.gene2).split(sep):
					if not gene_interval2:
						gene_interval2 = gene_ann.get_gene_interval(gene_name)
					else:
						break
				
			# Gene Strand : Defuse Strand
			# + >> + : + -
			# + >> - : + +
			# - >> + : - -
			# - >> - : - +
			if gene_interval1 and gene_interval2:
				if gene_interval1.strand == "f":
					self.strand1 = "+"
				else:
					self.strand1 = "-"
				if gene_interval2.strand == "f":
					self.strand2 = "-"
				else:
					self.strand2 = "+"
			else:
				# failed to fill in strand, return list with warnning info
				#gene_order.append("Strand_filling_failed")
				return gene_order

				
				
		if self.strand1 == "+" and self.strand2 == "+":
			gene_order = self.__check_gene_pairs(a, d, gene_ann)
			gene_order += self.__check_gene_pairs(b, c, gene_ann)
		elif self.strand1 == "+" and self.strand2 == "-":
			gene_order = self.__check_gene_pairs(a, b, gene_ann)
			gene_order += self.__check_gene_pairs(d, c, gene_ann)
		elif self.strand1 == "-" and self.strand2 == "+":
			gene_order = self.__check_gene_pairs(c, d, gene_ann)
			gene_order += self.__check_gene_pairs(b, a, gene_ann)
		elif self.strand1 == "-" and self.strand2 == "-":
			gene_order = self.__check_gene_pairs(c, b, gene_ann)
			gene_order += self.__check_gene_pairs(d, a, gene_ann)
		return gene_order

	# realign breakpoints of this fusion to the left most, not finished, how to define "left" when genes are on different chrs 
	def left_aln_fusion_bp(self, ref_file):
		rlen = 100
		refs = pysam.FastaFile(ref_file)

		#pysam use 0-based coordinates, cff use 1-based coordinates
		if self.strand1 == "+" and self.strand2 == "+":
			up_seq = sequtils.rc_seq(refs.fetch(self.chr1, self.pos1-rlen, self.pos1), "r")
			down_seq = sequtils.rc_seq(refs.fetch(self.chr2, self.pos2, self.pos2+rlen), "c")
		elif self.strand1 == "+" and  self.strand2 == "-":
			up_seq = sequtils.rc_seq(refs.fetch(self.chr1, self.pos1-rlen, self.pos1), "r")
			down_seq = sequtils.rc_seq(refs.fetch(self.chr2, self.pos2-rlen-1, self.pos2-1), "r")
		elif self.strand1 == "-" and  self.strand2 == "+":
			up_seq = refs.fetch(self.chr1, self.pos1-1, self.pos1+rlen-1)
			down_seq = refs.fetch(self.chr2, self.pos2, self.pos2+rlen)
		elif self.strand1 == "-" and  self.strand2 == "-":
			up_seq = sequtils.rc_seq(refs.fetch(self.chr1, self.pos1-1, self.pos1+rlen-1), "c")
			down_seq = sequtils.rc_seq(refs.fetch(self.chr2, self.pos2-rlen-1, self.pos2-1), "r")
		else:
			print >> sys.stderr, "Unknown strand:", self.strand1, self.strand2
			sys.exit(1)
		i = 0
		while i < rlen:
			if up_seq[i].upper() == down_seq[i].upper():
				i += 1
			else:
				break
		print up_seq.upper()
		print down_seq.upper()
		return i
				
class GeneIntervals():
	def __init__(self, bed_ann_list):
		'''
		self.gene_name = gene_name
		self.chr = chr
		self.start = start
		self.end = end
		self.strand = strand
		self.is_coding = is_coding
		'''
		self.load(bed_ann_list)
	def load(self, bed_ann_list):
		if not bed_ann_list:
			self.gene_name = "NA"
			self.chr = "NA"
			self.start = -1
			self.end = -1
			self.strand = "NA"
			self.is_coding = False
			#print >> sys.stderr, "Empty GeneBed annotation."
			#sys.exit(1)
		else:
			# warnning when same gene on different chr/strand, load  the first annotation
			if len(set([(a.gene_name, a.chr, a.strand) for a in bed_ann_list])) > 1:
				print >> sys.stderr, "Warnning: GeneBed annotation includes multiple genes."
				print >> sys.stderr, set([(a.gene_name, a.chr, a.strand) for a in bed_ann_list])
				#sys.exit(1)
			self.gene_name = bed_ann_list[0].gene_name
			self.chr = bed_ann_list[0].chr
			self.strand= bed_ann_list[0].strand

			self.start = min([a.start for a in bed_ann_list])
			self.end = max([a.end for a in bed_ann_list])
			self.is_coding = True if "cds" in [a.type for a in bed_ann_list] else False
				
	
	def overlap(self, interval2):
		if self.chr == interval2.chr and min(self.end, interval2.end) - max(self.start, interval2.end) > 0:
			return True
		else:
			return False
	def merge(self, interval2):
		if self.chr != interval2.chr or self.strand != interval2.strand:
			print >> sys.stderr, "Warning: intervals are on different chr/strand."
			print >> sys.stderr, self.gene_name, interval2.gene_name
			return self
		else:
			new_interval = GeneIntervals("Merged_" + self.gene_name + "_" + interval2.gene_name, self.chr, min(self,start, interval2.start), max(self.end, interval2.end), self.strand, False)
			return new_interval			
class BreakpointAnnotation:
	#tuple = (chr, start, end, transcript, type, idx, strand, gene_name)
	def __init__(self, chr, start, end, transcript_name, type, idx, strand, gene_name):
		self.__attrs = []
		self.chr = chr 			
		self.start = start
		self.end = end
		self.transcript_name = transcript_name
		self.type = type # utr/cds/intron
		self.idx = idx
		self.strand = strand
		self.gene_name = gene_name
		self.is_on_boundary = False
		self.__attrs.append(chr)
		self.__attrs.append(str(start))
		self.__attrs.append(str(end))
		self.__attrs.append(transcript_name)
		self.__attrs.append(type)
		self.__attrs.append(str(idx))
		self.__attrs.append(strand)
		self.__attrs.append(gene_name)
	def tostring(self):
		return  "\t".join(self.__attrs)
#Load bed format gene annotation, current support knowngene.bed's format, map given genomic loactions to genens, return matched gene list
class GeneAnnotation():
	__gene_starts = {}
	__genes = {}
	__max_diff = 1000000
	__gene_intervals = {}
	__transcripts_ann = {}
	__gene_name_id_map = {}
	__gene_name_idx_map = {}
	__coding_gene_list = []
	# gene_ann_bed
	# chr1    24683494	24685032	ENST00000003583 utr3    0       r       STPG1   ENSG00000001460
	def __init__(self, gene_ann_bed):
		if gene_ann_bed != "":
			self.load_gene_bed(gene_ann_bed)
			self.load_gene_intervals(gene_ann_bed)
			self.build_coding_gene_list()
			'''
			print >> sys.stderr, "mem of gene_starts:", sys.getsizeof(self.__gene_starts)
			print >> sys.stderr, "length:", len(self.__gene_starts.keys())
			print >> sys.stderr, "unit size:", sys.getsizeof(self.__gene_starts[self.__gene_starts.keys()[0]])
				
			print >> sys.stderr, "mem of genes:", sys.getsizeof(self.__genes)
			print >> sys.stderr, "length:", len(self.__genes.keys())
			print >> sys.stderr, "unit size:", sys.getsizeof(self.__genes[self.__genes.keys()[0]])

			print >> sys.stderr, "mem of gene_intervals:", sys.getsizeof(self.__gene_intervals)
			print >> sys.stderr, "length:", len(self.__gene_intervals.keys())
			print >> sys.stderr, "unit size:", sys.getsizeof(self.__gene_intervals[self.__gene_intervals.keys()[0]])
			'''
			#self.load_transcripts_ann(gene_ann_bed)
	# for each transcript save all its cds,intron, utr annotations in dictionary __transcripts_ann
	def load_transcripts_ann(self, gene_ann_bed):
		start_time = time.time()
		for line in open(gene_ann_bed, "r"):
			ann = GeneBed(line)
			key = ann.transcript_id
			self.__transcripts_ann.setdefault(key, []).append(ann)
		print >> sys.stderr, "Transcript annotations loaded."
		print >> sys.stderr, time.time() - start_time, "sec. elapsed."	
	
	def get_transcripts_ann(self, trans_id):
		if trans_id in self.__transcripts_ann:
			return self.__transcripts_ann[trans_id]
		else:
			return []	
	def get_gene_id(self, gene):
		if gene in self.__gene_name_id_map:
			return self.__gene_name_id_map[gene]
		else:
			return "NA"
	# for each gene use minimal start and max end of all its transcripts as its interval, same intervals of all genes in an dictionary __gnen_intervals
	def load_gene_intervals(self, gene_ann_bed):
		start_time = time.time()
		n1 = 0
		n2 = 0
		tmp_dict = {}
		# use gene name as key to bulid a dict for gene_ann_bed, all annotations for the same gene are saved in a list which will be used for getting intervals
		for line in open(gene_ann_bed, "r"):
			bed_ann = GeneBed(line)
			if self.filter_gene_ann(bed_ann):
				continue
			tmp_dict.setdefault(bed_ann.gene_name, []).append(bed_ann)
		for gene_name in tmp_dict:
			self.__gene_intervals.setdefault(gene_name, GeneIntervals(tmp_dict[gene_name]))
		print >> sys.stderr, "Gene intervals loaded."
		print >> sys.stderr, time.time() - start_time, "sec. elapsed."	
		#print >> sys.stderr, "#gene on diff chr:", n1, "#gene on diff strand:", n2

	# index foward/backward genes to decide whether two genes are next to each other (possible readthrough)
	def build_coding_gene_list(self):
		for gene_name in self.__gene_intervals:
			interval = self.__gene_intervals[gene_name]
			if interval.is_coding:
				self.__coding_gene_list.append(interval)
		self.__coding_gene_list.sort(key = lambda i:(i.chr, i.start))
		
		i_f = 0 # idx of forward starnd gene
		i_r = 0 # idx of reverse starnd gene
		#pre_interval_f = GeneIntervals("None", "chr0", 0, 0, "+", False)
		#pre_interval_r = GeneIntervals("None", "chr0", 0, 0, "-", False)
		pre_interval_f = GeneIntervals([])
		pre_interval_r = GeneIntervals([])
		for interval in self.__coding_gene_list:
			if interval.strand == "+" or interval.strand == "f":
				if interval.overlap(pre_interval_f):
					pre_interval_f = interval.merge(pre_interval_f)
				else:
					pre_interval_f = interval
					i_f += 1
				self.__gene_name_idx_map.setdefault(interval.gene_name, str(i_f) + "_" + interval.strand)
			elif interval.strand == "-" or interval.strand == "r":
				if interval.overlap(pre_interval_r):
					pre_interval_r = interval.merge(pre_interval_r)
				else:
					pre_interval_r = interval
					i_r += 1
				self.__gene_name_idx_map.setdefault(interval.gene_name, str(i_r) + "_" + interval.strand)

	# whether a gene include coding exon (cds)
	def is_coding(self, gene_name):
		return self.__gene_intervals[gene_name].is_coding
	
	
	def get_coding_gene_idx(self, gene_name):
		return self.__gene_name_idx_map[gene_name]

	def test_intervals(self):
		return self.__gene_intervals
	def filter_gene_ann(self, ann):
		if "hap" in ann.chr or "_" in ann.chr:
			return True
		else:
			return False

	def load_gene_bed(self, gene_ann_bed):
		start_time = time.time()
		n = 0
		for line in open(gene_ann_bed, "r"):
			ann = GeneBed(line)
			# filter out annotation not needed, e.g. gene annotation on chr_hap
			if self.filter_gene_ann(ann):
				continue
			self.__genes.setdefault(ann.chr, []).append(ann)
			# gene id and name map, this dict will take ~3G ram for ensgene annotation, current not loaded
			#self.__gene_name_id_map.setdefault(ann.gene_name, ann.gene_id)
			n += 1
		print >> sys.stderr, n, "annotations from", gene_ann_bed, "loaded."
		# sort all ann of each chr by start pos
		for chr in self.__genes:
			self.__genes[chr] = sorted(self.__genes[chr], key=lambda d:int(d.start))

			# save start pos in self.__gene_starts
			self.__gene_starts[chr] = [int(d.start) for d in self.__genes[chr]]
		print >> sys.stderr, time.time() - start_time, "sec. elapsed."	

	# for a given gene, return its interval if it is in the annotation, return empty interval if not.
	def get_gene_interval(self, gene_name):
		if gene_name in self.__gene_intervals:
			return self.__gene_intervals[gene_name]
		else:
			print >> sys.stderr, "Warnning: gene name", gene_name, "is not in current annotation."
			return None
	# get distance between two gene intervals, return -1 if genes are rom different chrs, minus value if overlap.
	def get_gene_distance(self, gene1, gene2):
		interval1 = self.get_gene_interval(gene1)
		interval2 = self.get_gene_interval(gene2)
		if interval1 and interval2:
			chr1 = interval1[0]
			start1 = int(interval1[1])
			end1 = int(interval1[2])
			chr2 = interval2[0]
			start2 = int(interval2[1])
			end2 = int(interval2[2])
			if chr1 != chr2: # genes from different chrs
				return -1 
			else:
				return cmp_intervals(start1, end1, start2, end2)
			
			
			
			
	def map_pos_to_genes(self, chr, pos):
		matched_genes = []
		if not chr in self.__gene_starts:
			return matched_genes
		idx = bisect.bisect(self.__gene_starts[chr], pos)	
		while 0 < idx <= len(self.__gene_starts[chr]):
			#bpann is an BreakpointAnnotation object
			bpann = self.__genes[chr][idx-1]
			#search within a limited region (default 1000000)
			if pos - bpann.start > self.__max_diff: 
				break

			if bpann.start <= pos <= bpann.end:
				if bpann.start == pos or bpann.end == pos:
					bpann.is_on_boundary = True
				matched_genes.append(bpann)
			idx -= 1
	
		return matched_genes
	def get_adjacent_exons(self, chr, pos):
		previous_exons = []
		next_exons = []	
			
		if not chr in self.__gene_starts:
			return previous_exons, next_exons
		idx = bisect.bisect(self.__gene_starts[chr], pos)
		while 0 < idx <= len(self.__gene_starts[chr]):
			
			bpann = self.__genes[chr][idx-1]
			#search within a limited region (default 1000000)
			if pos - bpann.start > self.__max_diff: 
				break
			if bpann.start <= pos <= bpann.end:
				# when breakpoint map to an intron/exon annotation
				#if bpann.type == "intron" or bpann.type == "cds": 
				if True: 
					# previous and next exons are based on coordinates, not strand. i.e. if a adjacent exon's coordinate < bp, it is a previous exon, otherwise next exon
					for i in range(idx, min(idx+100, len(self.__gene_starts[chr]))):
						adjacent_bpann = self.__genes[chr][i]
						# in next 100 annotations, try to find one exon next to current intron of the same transcript
						if bpann.transcript_name == adjacent_bpann.transcript_name and adjacent_bpann.type == "cds" and abs(adjacent_bpann.idx - bpann.idx) <= 3:
							next_exons.append(adjacent_bpann)
							#break
							#if bpann.strand == "f":
							#	next_exons.append(adjacent_bpann)
							#	break
							#else:
							#	previous_exons.append(adjacent_bpann)
							#	break
					for i in range(idx-2, max(0, idx-100), -1):
						adjacent_bpann = self.__genes[chr][i]
						# in next 100 annotations, try to find one exon next to current intron of the same transcript
						if bpann.transcript_name == adjacent_bpann.transcript_name and adjacent_bpann.type == "cds" and abs(adjacent_bpann.idx - bpann.idx) <= 3:
							previous_exons.append(adjacent_bpann)
							#break
							#if bpann.strand == "r":
							#	next_exons.append(adjacent_bpann)
							#	break
							#else:
							#	previous_exons.append(adjacent_bpann)
							#	break
			idx -= 1
		return previous_exons, next_exons	
		
		
			
				
class GeneFusions():
	__fusions = {}
	__genes = {}
	__geneann = GeneAnnotation("")
	def __init__(self, fusion_file, gene_ann_bed):
		self.__geneann = GeneAnnotation(gene_ann_bed)
		self.load_gene_fusions(fusion_file)
		
	# fusion_file is a tsv file in cff(common fusion format) format
	'''
	1 gene_chromosome1
	2 genomic_break_pos1
	3 genomic_strand1
	4 gene_chromosome2
	5 genomic_break_pos2
	6 genomic_strand2
	7 orf
	8 read_through
	9 splitr_count
	10 span_count
	11 Sample
	12 Run
	13 Tool
	14 id
	15 probability
	'''
	#The original genes reported by each tool are not included, will map the locations to gene to uniform the gene names
	def load_gene_fusions(self, fusion_file):
		for line in open(fusion_file, "r"):
			if line.startswith(("@", "#", "Chromosome1", "gene_")):
				continue

			fusion = CffFusion(line)
			gene_info1 = self.__geneann.map_pos_to_genes(fusion.chr1, fusion.pos1)
			gene_info2 = self.__geneann.map_pos_to_genes(fusion.chr2, fusion.pos2)
			gene_name1 = sorted(list(set([g[7] for g in gene_info1])))
			gene_name2 = sorted(list(set([g[7] for g in gene_info2])))
			gene_transcripts1 = ",".join(sorted(list(set([g[3] for g in gene_info1]))))
			gene_transcripts2 = ",".join(sorted(list(set([g[3] for g in gene_info2]))))

			# all mapped genes
			fusion.gene1 = ",".join(gene_name1)
			fusion.gene2 = ",".join(gene_name2)
			
			# all transcripts of each mapped gene
			fusion.trans_id1 = gene_transcripts1
			fusion.trans_id2 = gene_transcripts2

			# use gene pair of each fusion as key, for recurrent fusions
			if len(gene_name1)>0 and  len(gene_name2)>0:
				key = fusion.gene1 + "_" + fusion.gene2
				reverse_key = fusion.gene2 + "_" + fusion.gene1
				
				if not reverse_key in self.__fusions:
					self.__fusions.setdefault(key, []).append(fusion)
				else:
					self.__fusions.setdefault(reverse_key, []).append(fusion)

				# use one gene at a time as key, for recurrent genes
				for gene in fusion.gene1.split(",") + fusion.gene2.split(","):
					self.__genes.setdefault(gene, []).append(fusion)
			else:
				# some process for breakpoints mapped to one gene or no genes
				pass
	# get a set of sample names, remove duplicate samples		
	def get_samples_and_tools(self, fusion_list):
		#print [f.sample for f in fusion_list]
		return list(set([f.sample for f in fusion_list])), list(set([f.tool for f in fusion_list]))
	# get partner genes
	def get_partner_genes(self, gene, fusion_list):
		partner_genes = list(set([f.gene1 for f in fusion_list if f.gene1!=gene] + [f.gene2 for f in fusion_list if f.gene2!=gene]))
		return partner_genes
		
	def count_sample_number(self, sample_list):
		n_tumor = 0
		n_normal = 0
		for sample in sample_list:
			if sample.endswith("N"):
				n_normal += 1
			else:
				n_tumor += 1
		return n_normal, n_tumor
	def output_recurrent_fusions(self):
		# fusion looks like CHP1,MEP_TIMP2
		for fusion in self.__fusions:
			genes = fusion.split("_")
			genes1 = set(genes[0].split(","))
			genes2 = set(genes[1].split(","))
			mapped_to_same_gene = bool(genes1 & genes2)

			gene_ids1 = []
			gene_ids2 = []
			for gene in genes1:
				gene_ids1.append(self.__geneann.get_gene_id(gene))
			for gene in genes2:
				gene_ids2.append(self.__geneann.get_gene_id(gene))
		

			sample_fusion_list= self.__fusions[fusion]
			sample_list, tool_list = self.get_samples_and_tools(sample_fusion_list)
						
			n_normal, n_tumor = self.count_sample_number(sample_list)
			print "Gene_Fusion:", fusion, "Sample_num:", len(sample_list), "Gene_ID:", ",".join(gene_ids1) + "_" + ",".join(gene_ids2),  "Same_gene:", mapped_to_same_gene, "Tools:", "|".join(tool_list), "#Normal:", n_normal, "#Tumor", n_tumor, "Samples:", "|".join(sample_list)
			for sample_fusion in sample_fusion_list:
				print "\t" + sample_fusion.print_to_string()
				
	def output_recurrent_genes(self):
		for gene in self.__genes:
			sample_fusion_list= self.__genes[gene]
			sample_list, tool_list = self.get_samples_and_tools(sample_fusion_list)
			partner_genes = self.get_partner_genes(gene, sample_fusion_list)

			n_normal, n_tumor = self.count_sample_number(sample_list)
			
			
			print "Gene:", gene, "Sample_num:", len(sample_list), "#Normal:", n_normal, "#Tumor", n_tumor,  "Samples:", "|".join(sample_list), "Partner_num:", len(partner_genes), "Partner_genes:", "|".join(partner_genes)
			for sample_fusion in sample_fusion_list:
				print "\t" + sample_fusion.print_to_string()

	def test(self):
		self.output_recurrent_fusions()
		self.output_recurrent_genes()
		
		
			
#compare two genomic locations, return True if the difference is samller than parameter diff		
def cmp_breakpoints(chr1, pos1, chr2, pos2, diff):
	if chr1 != chr2:
		return False
	elif abs(pos1 - pos2) < diff:
		return True
# check whether two intervals overlap, overlap if end-start > 0
def cmp_intervals(start1, end1, start2, end2):
	start = max(start1, start2)
	end = min(end1, end2)
	return end - start


#return seq0's revers complement sequence
def rc_seq(seq0):
	seq=seq0.upper()
	tmpseq=seq.replace("A","a")
	seq=tmpseq
	tmpseq=seq.replace("T","A")
	seq=tmpseq
	tmpseq=seq.replace("a","T")
	seq=tmpseq
	tmpseq=seq.replace("G","g")
	seq=tmpseq
	tmpseq=seq.replace("C","G")
	seq=tmpseq
	tmpseq=seq.replace("g","C")
	seq=tmpseq[::-1]
	return seq

def pass_filter(read):
	if read.mapping_quality == 255: # tophat2 max mapq
		return True
	else:
		return False
# return #read per basepair for all exons(cds) of transcripts in trans_id_list
def get_exon_cov(trans_id_list, gene_ann, bam):

	# get uniquely mapped read number in each cds
	for trans_id in trans_id_list:
		for transcript in  gene_ann.get_transcripts_ann(trans_id):
			read_cnt = 0
			n = 0
			#if type == "cds":
			if 1:
				for read in bam.fetch(transcript.chr, transcript.start, transcript.end):
					n += 1
					# bam.fetch may return reads not in given region (bug of pysam?)
					# also require read to pass a "filter", currently only requires mapping quality to be maximum, 50 for tophat2, 60 for bwa, 255 for star
					if read.reference_start in range(transcript.start, transcript.end) and pass_filter(read):
						read_cnt += 1
					else:
						pass
						#print read.reference_start, read.reference_start in [start, end]
			
				print transcript.tostring, read_cnt, float(read_cnt)/(transcript.end - transcript.start)
def get_fusion_exon_cov(gene_fusion, gene_ann, bam):
	# map fusion breakpoints to genes
	matched_genes1 = gene_ann.map_pos_to_genes(gene_fusion.chr1, gene_fusion.pos1)	
	matched_genes2 = gene_ann.map_pos_to_genes(gene_fusion.chr2, gene_fusion.pos2)	
	
	# get all transcripts of genes
	trans_id_list1 = list(set(g[3] for g in matched_genes1))
	trans_id_list2 = list(set(g[3] for g in matched_genes2))
	trans_id_list = trans_id_list1 + trans_id_list2
	get_exon_cov(trans_id_list, gene_ann, bam)

#build fusion reference, refs is pysam.FastaFile, seg_len is the length of extraced ref sequence
def output_fusion_fasta(fusion, refs, seg_len):
	## defuse strands have different meanings, need to convert to the following pattern before uesd
	# ++ : ->| |->
	# +- : ->| <-|
	# -+ : |<- |->
	# -- : |<- <-|
	## defuse:
	# ++ : ->| <-|
	# +- : ->| |->
	# -+ : |<- <-|
	# -- : |<- |->
	
	# pysam.FastaFile uses 0-based coordinates, cff fusion uses 1-based
	# strand used here is differnt from defuse
	if strand1 == "+": ## ->|
		win_start1 = fusion.pos1 - seg_len
		win_end1 = fusion.pos1
		seq1 = refs.fetch(fusion.chr1, win_start1, win_end1)
	else: ## |<-
		win_start1 = fusion.pos1 - 1
		win_end1 = fusion.pos1 + seg_len - 1
		seq1 = rc_seq(refs.fetch(fusion.chr1, win_start1, win_end1))

	if strand2 == "+": ## |->
		win_start2 = fusion.pos2 - 1
		win_end2 = fusion.pos2 + seg_len - 1
		seq2 = refs.fetch(fusion.chr2, win_start2, win_end2)
	else: ## <-|
		win_start2 = fusion.pos2 - seg_len
		win_end2 = fusion.pos2
		seq2 = rc_seq(refs.fetch(fusion.chr2, win_start2, win_end2))
	seq = seq1 + seq2
# build transcript sequence based on bed format annotation file
# chr1    24683494	24685032	ENST00000003583 utr3    0       r       STPG1   ENSG00000001460
	
def build_transcript_seq(gene_ann, ref):
	trans_seqs = {}
	for chr in gene_ann.__transcripts_ann:
		seq = []
		for transcript in gene_ann.__transcripts_ann[chr]: 	
			seq.append(refs.fetch(transcript.chr, transcript.start, transcript.end))		
				
		trans_seqs.setdefault(transcript, "".join(seq))
	return trans_seqs
### old code	
# load bed format gene annotations as an dict based index
def load_gene_bed(bed_file):
	
	dict_starts = {}
	dict_genes = {}
	for line in open(bed_file, "r"):
		tmp = line.split()
		chr = tmp[0]
		start = int(tmp[1]) + 1 # bed coordinate is 0-based
		end = int(tmp[2])
		transcript = tmp[3]
		type = tmp[4]
		id = int(tmp[5])
		strand = tmp[6]
		gene = tmp[7]


		tuple = (start, end, transcript, type, id, strand, gene)
		dict_genes.setdefault(chr, []).append(tuple)
	# sort all ann of each chr by start pos
	for chr in dict_genes:
		dict_genes[chr] = sorted(dict_genes[chr], key=lambda d:int(d[0]))
		# save start pos in dict_starts
		dict_starts[chr] = [int(tuple[0]) for tuple in dict_genes[chr]]
	return dict_starts, dict_genes

# check if there are indels/clippings in cigar
def has_indels_or_clippings(cigar):
	for s in ["S","H","I","D","N"]:
		if s in cigar:
			return True
	return False

def print_read(read, type, outf):
	if not type in ["sam", "fq", "fa"]:
		print >> sys.stderr, "Unknown output type:", type
		sys.exit(1)
	name = read.query_name
	seq = read.query_sequence
	qual = "I" * len(seq) ## pysam provides number format qualities, to be fixed
	if read.is_reverse:
		seq = rc_seq(seq)

	if type == "sam":
		print read.tostring(samfile),
	elif type == "fq":
		print >> outf, "@"+name
		print >> outf, seq
		print >> outf, "+"
		print >> outf, qual
	elif type == "fa":
		print >> outf, ">" + name
		print >> outf, seq
	#print read.tostring(samfile),


# each read_list contains all reads (i.e. all alingments of a read pair) with the same name, descide whether to output this pair
def map_pos_to_genes_bak(read):
	matched_genes = []
	chr = "chr" + read.reference_name
	if not chr in dict_starts:
		return matched_genes
	pos = read.reference_start + 1 # reference_start is 0-based
	idx = bisect.bisect(dict_starts[chr], pos)	
	while 0 < idx < len(dict_starts[chr]):
		#tuple structure: tuple = (start, end, transcript, type, id, strand, gene)
		#tuple = dict_genes[chr][idx-1]
		#start = tuple[0]
		bpann = dict_genes[chr][idx-1]
		start = bpann.start

		#check within a limited region (100000)
		if pos - start > 100000: 
			break

		if bpann.start <= pos <= bpann.end:
			matched_genes.append(bpann)
			idx -= 1
		else:
			break
	return matched_genes

def is_concordant(read1, read2):
	if read1.reference_name == read2.reference_name and \
		((read1.reference_start < read2.reference_start and not read1.is_reverse and read2.is_reverse) or \
		(read2.reference_start <= read1.reference_start and not read2.is_reverse and read1.is_reverse)):
		return True
	else:
		return False
def filter_read(read_list, dict_starts, dict_genes, outf_1, outf_2): 
	if len(read_list) == 0:
		return

	read1 = ""
	read2 = ""

	f_unmap = False	# either read unmapped
	f_mapq = False	# either mapq == 60
	f_cigar = False	# no indel/clippings in cigar
	f_mismatch = False # miscmatch <= 1 
	f_gene = False	# both reads mapped to exon/utr
	f_concordant = False # check strand and order, not insert size
	n = 0
	for read in read_list:
		if read.is_secondary or read.is_supplementary:
			continue
		else:
			if read.is_read1:
				n += 1
				read1 = read
			else:
				n += 10
				read2 = read
	if n != 11:
		#print >> sys.stderr, "Read name", read1.query_name, "contians", n, "reads."
		print >> sys.stderr, "err1", str(read_list)
		sys.exit(1)
#	if read1 == "" or len(read2) == "":
		#print >> sys.stderr, "Read name ", read1.query_name, " is not paired."
#		print >> sys.stderr, "err2", str(read_list)
#		sys.exit(1)
	#either or both reads unaligned
	if read1.is_unmapped or read2.is_unmapped:
		f_unmap = True
	else:
		# either read contains >1 mismatch
		if read1.get_tag("NM") > 1 or read2.get_tag("NM") > 1:
			f_mismatch = True
		# both reads uniquely aligned
		if read1.mapping_quality == 60 and read2.mapping_quality == 60:
			f_mapq = True
		# both reads mapped continuously
		if (has_indels_or_clippings(read1.cigarstring) or has_indels_or_clippings(read2.cigarstring)):
			f_cigar = True
		
		# read pair mapped discordantly (strand and order)
		
		if not is_concordant(read1, read2):
			f_concordant = True

		# either read not mapped to exons/utrs
		matched_genes1 = set([g[6] for g in map_pos_to_genes(read1)])
		matched_genes2 = set([g[6] for g in map_pos_to_genes(read2)])
		
		if len(matched_genes1) == 0 or len(matched_genes2) == 0:
			f_gene = True
		elif len(matched_genes1.intersection(matched_genes2)) == 0:
			#print "set1", matched_genes1
			#print "set2", matched_genes2
			#print_read(read1, "sam", "")
			#print_read(read2, "sam", "")
			f_gene = True
	'''		
	print read1.cigarstring, read2.cigarstring
	print "f_unmap:", f_unmap
	print "f_mapq:", f_mapq
	print "f_cigar:", f_cigar
	print "f_gene:", f_gene
	print "f_mismatch:", f_mismatch
	'''
	if f_unmap or (f_mapq and (f_cigar or f_gene or f_mismatch or f_concordant)):
		print_read(read1, "fq", outf_1)
		print_read(read2, "fq", outf_2)

	#q = pybedtools.create_interval_from_list(["chr" + read1.reference_name, str(read.query_alignment_start), str(read.query_alignment_start)])
	#print(exon_bed_intervals.any_hits(q))
	#print(exon_bed_intervals.all_hits(q))
	
#print  >> sys.stderr, "load_gene_model"
