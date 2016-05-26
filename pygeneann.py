#!/usr/bin/env python

import sys
import pysam
import bisect
import time
import sequtils
import copy
print >> sys.stderr, "Version 0.2"
class GeneBed():
	def __init__(self, bed_line):
		tmp = bed_line.split()
		self.bed_line = bed_line.strip()
		self.chr = tmp[0]
		self.start = int(tmp[1]) + 1	# bed file use 0-based coordinate
		self.end = int(tmp[2])		# start and end are first and last base of each segment
		self.transcript_id = tmp[3]
		self.type = tmp[4] # utr/intron/cds
		self.idx = int(tmp[5])
		self.strand = tmp[6]
		self.gene_name = tmp[7]
		self.gene_id = tmp[8]
		# the following two attr are only used for breakpoint annotation, tell whether current breakpoint on/close to boundary, need to reset everytime before used
		self.is_on_boudary = False
		self.close_to_boundary = False
		# in gene intervals, this gene includes annotations from different chr, strand or >5Mb distance
		self.is_from_contradictory_gene = False
	def tostring(self):
		return self.bed_line
class CffFusionStats():
	__fusion_dict = {}
	__fusion_samples_dict = {}
	def __init__(self, cff_file):
		pass
		#self.__load_fusions(cff_file)
	def __load_fusions(self, cff_file):
		for line in open(cff_file, "r"):
			fusion = CffFusion(line)
			key = fusion.t_gene1 + "_" + fusion.t_gene2
			rkey = fusion.t_gene2 + "_" + fusion.t_gene1
			if rkey not in self.__fusion_dict:
				self.__fusion_dict.setdefault(key, []).append(fusion)
				self.__fusion_samples_dict.setdefault(key, []).append(fusion.sample_name)
			else:
				self.__fusion_dict.setdefault(rkey,[]).append(fusion)
				self.__fusion_samples_dict.setdefault(rkey, []).append(fusion.sample_name)
	
	# compare two fusions based on the up/downstream gene sets, if overlap on both sets consider them as same fusion
	def is_same_gene_pair_fusion(self, fusion1, fusion2):
		id = 1 
		# compare gene sets on fw strand 
		up_g1, down_g1 = set(fusion1.get_reannotated_genes(id))
		up_g2, down_g2 = set(fusion2.get_reannotated_genes(id))

		if up_g1 & up_g2 and down_g1 & down_g2:
			return True

		id = 2 
		# compare gene sets on bw strand 
		up_g1, down_g1 = set(fusion1.get_reannotated_genes(id))
		up_g2, down_g2 = set(fusion2.get_reannotated_genes(id))

		if up_g1 & up_g2 and down_g1 & down_g2:
			return True

		return False
	# compare gene fusions based on breakpoints
	def generate_common_fusion_stats(self, cff_file):
		fusion_bp_dict = {}
		fusion_pos1_idx_dict = {}
		# build index for fusions based on fusion.pos1
		for line in  open(cff_file, "r"):
			if line.startswith("#"):
				continue
			fusion = CffFusion(line)
			fusion_bp_dict.setdefault(fusion.chr1, {}).setdefault(fusion.strand1, []).append(fusion)
		for chr in fusion_bp_dict:
			for strand in fusion_bp_dict[chr]:
				fusion_bp_dict[chr][strand].sort(key = lambda x:x.pos1)	
				fusion_pos1_idx_dict.setdefault(chr, {}).setdefault(strand, [f.pos1 for f in fusion_bp_dict[chr][strand]]) 
		# find common fusions
	# output fusions in a fusion list as clustered fusions
	def output_clustered_fusions(self, fusion_list, cluster_type):	
			sample_list = [f.sample_name for f in fusion_list]	
			disease_list = [f.disease for f in fusion_list]	
			tool_list = [f.tool for f in fusion_list]	
			sample_type_list = [f.sample_type for f in fusion_list]	

			category_list = [f.reann_category1 for f in fusion_list]	
			category_list += [f.reann_category2 for f in fusion_list]	
			gene_order_list = [f.reann_gene_order1 for f in fusion_list]
			gene_order_list += [f.reann_gene_order2 for f in fusion_list]
			print cluster_type, ",".join(list(set(sample_list))), ",".join(list(set(sample_type_list))), ",".join(list(set(disease_list))), ",".join(list(set(tool_list))), ",".join(list(set(category_list))), ";".join(list(set(gene_order_list)))
	# cluster fusions of "NoDriverGene" and "Truncated" type on their breakpoints
	def generate_common_fusion_stats_by_breakpoints(self, fusion_list):
		diff = 100000
		# save clustered fusion id, skip a fusion when it is already clustered
		clustered_id = {}
		print >> sys.stderr, "bp fusion list:", len(fusion_list)
		for i in range(len(fusion_list)):
			if i in clustered_id:
				continue

			#if i % 100 == 0: print >> sys.stderr, i

			fusion1 = fusion_list[i]
			small_bp1, big_bp1 = fusion1.get_ordered_breakpoints()
			clustered_id.setdefault(i, i)
			fusion_cluster_list = []
			fusion_cluster_list.append(fusion1)
			for j in range(len(fusion_list)):
				if j in clustered_id:
					continue
				if i == j:
					continue

				fusion2 = fusion_list[j]
				small_bp2, big_bp2 = fusion2.get_ordered_breakpoints()
				if cmp_fusion_breakpoints(small_bp1, small_bp2, diff) and cmp_fusion_breakpoints(big_bp1, big_bp2, diff):
					clustered_id.setdefault(j, j)
					fusion_cluster_list.append(fusion2)

			self.output_clustered_fusions(fusion_cluster_list, "BP_Cluster")
				

	# compare gene fusions based on re-annotated gene names, hard to work for truncated fusions, used generate_common_fusion_stats_by_breakpoints for NoDriverGenes and Trucated type fusions
	def generate_common_fusion_stats_by_genes(self, cff_file):
		fusion_dict = {}
		fusion_list_for_bp_cmp = []
		common_key_dict = {}
		for line in  open(cff_file, "r"):
			if line.startswith("#"):
				continue
			fusion = CffFusion(line)
			keys = set()
			found = False # found both up/down stream genes
			id = 1
			if fusion.reann_category1 not in ["SameGene", "NoDriverGene"]:
				up_g, down_g = fusion.get_reannotated_genes(id)
				for g1 in up_g:
					for g2 in down_g:
						key = g1 + "_" + g2 
						keys.add(key)
						fusion_dict.setdefault(key, []).append(fusion)
						found = True
			id = 2
			if fusion.reann_category2 not in ["SameGene", "NoDriverGene"]:
				up_g, down_g = fusion.get_reannotated_genes(id)
				for g1 in up_g:
					for g2 in down_g:
						key = g1 + "_" + g2 
						keys.add(key)
						fusion_dict.setdefault(key, []).append(fusion)
						found = True
			
			if not found:
				fusion_list_for_bp_cmp.append(fusion)
			# for breakpoints that can be mapped to multiple genes pairs, save these pairs in a dict, merge their fusion lists when output
			if len(keys) > 1:
				for key in keys:
					if key in common_key_dict:
						common_key_dict[key] |= keys
					else:
						common_key_dict[key] = keys
		removed_keys = {} 
		for key in fusion_dict:
			if key in removed_keys:
				continue
			fusion_list = fusion_dict[key]

			if key in common_key_dict:
				key_set = common_key_dict[key]
				while True:
					l = len(key_set)
					key_set2 = key_set.copy()
					for key2 in key_set2:
						if key2 != key:
							key_set |= common_key_dict[key2]	
					if l == len(key_set):
						break
				#print key, key_set
				key_set.remove(key)
				for key2 in key_set:
					fusion_list += fusion_dict[key2]
					removed_keys.setdefault(key2, key2)	
			self.output_clustered_fusions(fusion_list, "Gene_Cluster")
			'''
			sample_list = [f.sample_name for f in fusion_list]	
			disease_list = [f.disease for f in fusion_list]	
			tool_list = [f.tool for f in fusion_list]	
			sample_type_list = [f.sample_type for f in fusion_list]	

			category_list = [f.reann_category1 for f in fusion_list]	
			category_list += [f.reann_category2 for f in fusion_list]	
			gene_order_list = [f.reann_gene_order1 for f in fusion_list]
			gene_order_list += [f.reann_gene_order2 for f in fusion_list]
			print key, ",".join(list(set(sample_list))), ",".join(list(set(sample_type_list))), ",".join(list(set(disease_list))), ",".join(list(set(tool_list))), ",".join(list(set(category_list))), ",".join(list(set(gene_order_list)))
			'''
		# send "NoDriverGene" and "Truncated" type fusions for breakpoint cluster
		self.generate_common_fusion_stats_by_breakpoints(fusion_list_for_bp_cmp)		
	def get_gene_order_stats(self):
		for key in self.__fusion_dict:
			n_sg = 0 # same gene
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
				#tmp = ";".join(fusion.otherann)
				tmp = []
				for attr in fusion.zone4_attrs:
					tmp.append(fusion.__dict__[attr])
				gene_order.append(" ".join(tmp))
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
				print "Fusion", fusion.t_gene1, fusion.t_gene2, ",".join(list(set(self.__fusion_samples_dict[key]))), type, ",".join(sorted(list(set(sample_type)))),
				print "|".join(list(set(gene_order)))
				print "\tSameGene:", n_sg
				print "\tReadThrough:", n_rt
				print "\tGeneFusion:", n_gf
				print "\tTruncatedCoding:", n_tc
				print "\tTruncatedNoncoding:", n_tn
				print "\tNonSense:", n_ns


class CffFusion():
	# chr1 pos1 strand1 chr2 pos2 strand2 library sample_name sample_type disease tool split_cnt span_cnt tool_gene1 tool_gene_area1 tool_gene2 tool_gene_area2 fwd_gene_order bwd_gene_order
	def __init__(self, cff_line):
		tmp = cff_line.split()
		# Breadkpoint Zone
		self.chr1 = tmp[0]
		self.pos1 = int(tmp[1])
		self.strand1 = tmp[2]
		self.chr2 = tmp[3]
		self.pos2 = int(tmp[4])
		self.strand2 = tmp[5]
		# Sample info Zone
		self.library = tmp[6] # DNA/RNA
		self.sample_name = tmp[7]
		self.sample_type = tmp[8] # Tumor/Normal
		self.disease = tmp[9]
		# Software Zone
		self.tool = tmp[10]
		self.split_cnt = int(tmp[11])
		self.span_cnt = int(tmp[12]) if tmp[12] != "NA" else tmp[12]
		self.t_gene1 = tmp[13] # gene reported by tool
		self.t_area1 = tmp[14] # exon/utr/intron
		self.t_gene2 = tmp[15]
		self.t_area2 = tmp[16]
		# Re-annotation Zone
		if len(tmp) >= 25:
			self.reann_gene_order1 = tmp[17] # re-annotated gene order on fw strand e.g. ZNF248_utr3,cds>>RP11-258F22.1_utr3
			self.reann_gene_type1 = tmp[18] # re-annotated gene type e.g. CodingGene>>NoncodingGene
			self.reann_gene_index1 = tmp[19] # gene index for read-through inference e.g. 7409_r>>9974_r, only valid for coding gene
			self.reann_category1 = tmp[20] # infered fusion type, ReadThrough, GeneFusion, TruncatedCoding, TruncatedNoncoding, Nosenes, SameGene

			self.reann_gene_order2 = tmp[21] # backward re-annotation
			self.reann_gene_type2 = tmp[22] 
			self.reann_gene_index2 = tmp[23] 
			self.reann_category2 = tmp[24]

		else:
			self.reann_gene_order1 = "NA" # re-annotated gene order on fw strand e.g. ZNF248_utr3,cds>>RP11-258F22.1_utr3
			self.reann_gene_type1 = "NA" # re-annotated gene type e.g. CodingGene>>NoncodingGene
			self.reann_gene_index1 = "NA" # gene index for read-through inference e.g. 7409_r>>9974_r, only valid for coding gene
			self.reann_category1 = "NA" # infered fusion type, ReadThrough, GeneFusion, TruncatedCoding, TruncatedNoncoding, Nosenes, SameGene

			self.reann_gene_order2 = "NA" # backward re-annotation
			self.reann_gene_type2 = "NA"
			self.reann_gene_index2 = "NA"
			self.reann_category2 = "NA"
		self.boundary_info = ""		
		# same all attrs in a list, for printing	
		self.zone1_attrs = ["chr1", "pos1", "strand1", "chr2", "pos2", "strand2"]
		self.zone2_attrs = ["library", "sample_name", "sample_type", "disease"]
		self.zone3_attrs = ["tool", "split_cnt", "span_cnt", "t_gene1", "t_area1", "t_gene2", "t_area2"]
		self.zone4_attrs = ["reann_gene_order1", "reann_gene_type1", "reann_gene_index1", "reann_category1", "reann_gene_order2", "reann_gene_type2", "reann_gene_index2", "reann_category2"]
		# format chr
		if not self.chr1.startswith("chr"):
			self.chr1 = "chr" + self.chr1
		if not self.chr2.startswith("chr"):
			self.chr2 = "chr" + self.chr2
		# check fields
		if not self.check_cff():
			sys.exit(1)
	def get_gene_names_from_gene_order(self):
		g1 = []
		g2 = []
		g3 = []
		g4 = []
		if self.reann_gene_order1 != "NA":
			tmp = self.reann_gene_order1.split(">>")
			tmp2 = tmp[0].split(";")
			for t2 in tmp2:
				tmp3 = t2.split("_")
				g1.append(tmp3[0])
			tmp2 = tmp[1].split(";")
			for t2 in tmp2:
				tmp3 = t2.split("_")
				g2.append(tmp3[0])
			
			
		if self.reann_gene_order2 != "NA":
			tmp = self.reann_gene_order2.split(">>")
			tmp2 = tmp[0].split(";")
			for t2 in tmp2:
				tmp3 = t2.split("_")
				g3.append(tmp3[0])
			tmp2 = tmp[1].split(";")
			for t2 in tmp2:
				tmp3 = t2.split("_")
				g4.append(tmp3[0])
		
		return g1, g2, g3, g4		

	def __test_prefered_gene(self, genes, gene_name_list):
		prefered_gene = ""
		for gene in genes:
			if gene.gene_name in gene_name_list:
				if not prefered_gene:
					prefered_gene = gene
				else:
					if gene.type == "cds":
						prefered_gene = gene
						break
					elif gene.type == "intron" and prefered_gene != "cds":
						prefered_gene = gene
						
		return prefered_gene		
			
	# among all the re-annotated genes, return a prefered one which should be: 1. breakpoint in/on exon. 2. if not, has a exon closest to breakpoint. 3. the longest isoform.
	# Not finished, needs to consider
	def get_perfered_ann(self, gene_ann):
		print >> sys.stderr,  "Not finished."
		return
		genes1 = gene_ann.map_pos_to_genes(self.chr1, self.pos1)
		genes2 = gene_ann.map_pos_to_genes(self.chr2, self.pos2)
		# get gene name list from re-annotated gene orders
		g1, g2, g3, g4 = self.get_gene_names_from_gene_order()
		
		prefered_gene1 = self.__test_prefered_gene(genes1, g1)
		prefered_gene2 = self.__test_prefered_gene(genes2, g2)
		# try the opposite, because breakpoints order are not always consistant with gene orders
		if not (prefered_gene1 and prefered_gene2):
			prefered_gene1 = self.__test_prefered_gene(genes1, g2)
			prefered_gene2 = self.__test_prefered_gene(genes2, g1)
			
		prefered_gene3 = self.__test_prefered_gene(genes1, g3)
		prefered_gene4 = self.__test_prefered_gene(genes2, g4)
			
		if not (prefered_gene3 and prefered_gene4):
			prefered_gene3 = self.__test_prefered_gene(genes1, g4)
			prefered_gene4 = self.__test_prefered_gene(genes2, g3)
		
		print self.tostring()
		if prefered_gene1 and prefered_gene2:
			print "\t" + prefered_gene1.tostring()
			print "\t" + prefered_gene2.tostring()
		elif prefered_gene3 and prefered_gene4:		
			print "\t" + prefered_gene3.tostring()
			print "\t" + prefered_gene4.tostring()
		else:
			print "No prefered genes."
		
	# compare fusion breakpoints and return in an order of smaller to  bigger
	def get_ordered_breakpoints(self):
		if self.chr1 < self.chr1:
			small_bp = (self.chr1, self.pos1, self.strand1)
			big_bp = (self.chr2, self.pos2, self.strand2)
		elif self.chr1 == self.chr1 and self.pos1 < self.pos2:
			small_bp = (self.chr1, self.pos1, self.strand1)
			big_bp = (self.chr2, self.pos2, self.strand2)
		else:
			small_bp = (self.chr2, self.pos2, self.strand2)
			big_bp = (self.chr1, self.pos1, self.strand1)
		return small_bp, big_bp	
				
	# after reannotation, get a list of upstream genes and a list of downstream genes, can be used to compare fusions from different tools on gene leverl
	def get_reannotated_genes(self, id):
		up_genes = []
		down_genes = []
		if id == 1:
			gene_order = self.reann_gene_order1
		elif id == 2:
			gene_order = self.reann_gene_order2
		else:
			print >> sys.stderr, "Strand has to be 1 or 2", id, "provided."
			sys.exit(1)
			
		#LINC00875_utr5,intron>>NBPF9_intron;LOC100288142_intron;NBPF8_intron    NoncodingGene>>CodingGene,CodingGene,CodingGene >>566_f,525_r,565_f     TruncatedNoncoding
		if gene_order != "NA":
			tmp = gene_order.split(">>")
			for g in tmp[0].split(";"):
				gname = g.split("_")[0]
				if g: up_genes.append(gname)
			for g in tmp[1].split(";"):
				gname = g.split("_")[0]
				if g: down_genes.append(gname)
		
		if not up_genes:
			up_genes.append("InterGenic")		
		if not down_genes:
			down_genes.append("InterGenic")		
		return set(up_genes), set(down_genes)
	def check_cff(self):
		if self.library not in ["NA", "DNA", "RNA"]:
			print >> sys.stderr, "Unknown library type:", self.library
			return False
		if self.sample_type not in ["NA", "Tumor", "Normal"]:
			print >> sys.stderr, "Unknown sample type:", self.sample_type
			return False
		return True

	def tostring(self):
		value = []
		for attr in self.zone1_attrs + self.zone2_attrs + self.zone3_attrs + self.zone4_attrs:
			if not attr in self.__dict__:
				print >> sys.stderr, "Attribute name error:", attr
				sys.exit(1)
			else:
				value.append(self.__dict__[attr])
		return "\t".join(map(lambda x:str(x), value)) + "\t" + self.boundary_info
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

		list1 = []
		for gene_name in genes1:
			tmp = set(genes1[gene_name])
			if gene_ann.is_contradictory(gene_name):
				gene_name += "(Cont)"
			list1.append(gene_name + "_" + ",".join(list(tmp)))
		list2 = []
		for gene_name in genes2:
			tmp = set(genes2[gene_name])
			if gene_ann.is_contradictory(gene_name):
				gene_name += "(Cont)"
			list2.append(gene_name + "_" + ",".join(list(tmp)))
		# No driver gene
		if not genes1:
			#gene_order.extend(["NA"]*3)
			#gene_order.append("NoDriverGene")
			category = "NoDriverGene"
		# map to same gene
		elif common_genes:
			#gene_order.append("SameGene" + "\t" + ",".join(list(common_genes)))
			#gene_order.extend(["NA"]*3)
			#gene_order.append("SameGene")
			category = "SameGene"
		else:
			#gene_order.append(",".join(list(set(genes1))) + ">>" + ",".join(list(set(genes2))))

			# category fusions into: read through, gene fusion, truncated coding, truncated noncoding, nonsense
			if "CodingGene" in type1 and "CodingGene" in type2:
				for id in id1:
					tmp = id.split("_")
					idx = int(tmp[0])
					strand = tmp[1]
					#ReadThrough: gene1 and gene2 are adjacent genes (id1 - id2 = 1) or overlapping genes (id1 = id2) but breakpoints cannot map to same gene
					if (strand == "f" and  str(idx+1) + "_f" in id2 ) or (strand == "r" and  str(idx-1) + "_r" in id2) or (id in id2):
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



		gene_order.append(",".join(list1) + ">>" + ";".join(list2))
		gene_order.append(",".join(type1) + ">>" + ",".join(type2))
		gene_order.append(",".join(id1) + ">>" + ",".join(id2))
		gene_order.append(category)

		return gene_order
	# based on given gene annotations re-annotate cff fusions, infer possible up/downstream genens, try to fill in strand if info missing
	def ann_gene_order(self, gene_ann):
		gene_order = []
		# fusion has been annotated
		if self.reann_category1 != "NA" or self.reann_category2 !="NA":
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
				for gene_name in (self.t_gene1).split(sep):
					if not gene_interval1:
						gene_interval1 = gene_ann.get_gene_interval(gene_name)
					else:
						break
				for gene_name in (self.t_gene2).split(sep):
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

				
		# gene_order includes: 5' gene >> 3' gene, 5' gene type >> 3' gene type, 5' coding gene idx >> 3' coding gene inx, category
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
		self.reann_gene_order1, self.reann_gene_type1, self.reann_gene_index1, self.reann_category1, self.reann_gene_order2, self.reann_gene_type2, self.reann_gene_index2, self.reann_category2 = gene_order
		
		# check whether pos on mapped genes boundaries
		on_boundary1 = False
		close_to_boundary1 = False
		on_boundary2 = False
		close_to_boundary2 = False
		for gene in matched_genes1:
			if not gene.type == "cds":
				break
			if gene.is_on_boundary:
				on_boundary1 = True
				close_to_boundary1 = True
				break
			if gene.close_to_boundary:
				close_to_boundary1 = True
		for gene in matched_genes2:
			if not gene.type == "cds":
				break
			if gene.is_on_boundary:
				on_boundary2 = True
				close_to_boundary2 = True
				break
			if gene.close_to_boundary:
				close_to_boundary2 = True
		self.boundary_info = "\t".join(map(str, [on_boundary1, close_to_boundary1, on_boundary2, close_to_boundary2]))
		gene_order += [on_boundary1, close_to_boundary1, on_boundary2, close_to_boundary2]
		return gene_order

	# realign breakpoints of this fusion to the left most, not finished, how to define "left" when genes are on different chrs 
	def left_aln_fusion_bp(self, refs):
		# provided reference file lacks fusion chr
		if not (self.chr1 in refs.references and self.chr2 in refs.references):
			return (-1, -1)
		rlen = 10
		#refs = pysam.FastaFile(ref_file)
		
		#pysam use 0-based coordinates, cff use 1-based coordinates
		if self.strand1 == "+" and self.strand2 == "+":
			#up_seq = sequtils.rc_seq(refs.fetch(self.chr1, self.pos1-rlen, self.pos1+rlen), "r")
			#down_seq = sequtils.rc_seq(refs.fetch(self.chr2, self.pos2-rlen, self.pos2+rlen), "c")
			up_seq = refs.fetch(self.chr1, self.pos1-rlen, self.pos1+rlen)
			down_seq = sequtils.rc_seq(refs.fetch(self.chr2, self.pos2-rlen, self.pos2+rlen), "rc")
		elif self.strand1 == "+" and  self.strand2 == "-":
			#up_seq = sequtils.rc_seq(refs.fetch(self.chr1, self.pos1-rlen, self.pos1+rlen), "r")
			#down_seq = sequtils.rc_seq(refs.fetch(self.chr2, self.pos2-rlen-1, self.pos2-1+rlen), "r")
			up_seq = refs.fetch(self.chr1, self.pos1-rlen, self.pos1+rlen)
			down_seq = refs.fetch(self.chr2, self.pos2-rlen-1, self.pos2-1+rlen)
		elif self.strand1 == "-" and  self.strand2 == "+":
			#up_seq = refs.fetch(self.chr1, self.pos1-1-rlen, self.pos1+rlen-1)
			#down_seq = refs.fetch(self.chr2, self.pos2-rlen, self.pos2+rlen)
			down_seq = refs.fetch(self.chr1, self.pos1-1-rlen, self.pos1+rlen-1)
			up_seq = refs.fetch(self.chr2, self.pos2-rlen, self.pos2+rlen)
		elif self.strand1 == "-" and  self.strand2 == "-":
			#up_seq = sequtils.rc_seq(refs.fetch(self.chr1, self.pos1-1-rlen, self.pos1+rlen-1), "c")
			#down_seq = sequtils.rc_seq(refs.fetch(self.chr2, self.pos2-rlen-1, self.pos2-1+rlen), "r")
			#up_seq = sequtils.rc_seq(refs.fetch(self.chr1, self.pos1-rlen, self.pos1+rlen), "rc")
			#down_seq = refs.fetch(self.chr2, self.pos2-rlen-1, self.pos2-1+rlen)
			down_seq = refs.fetch(self.chr1, self.pos1-rlen-1, self.pos1+rlen-1)
			up_seq = sequtils.rc_seq(refs.fetch(self.chr2, self.pos2-rlen-1, self.pos2+rlen-1), "rc")
		else:
			print >> sys.stderr, "Unknown strand:", self.strand1, self.strand2
			sys.exit(1)
		
		if len(up_seq) < rlen or len(down_seq) < rlen:
			print >> sys.stderr, "Warnning: reference sequence cannot be fetched."
			print >> sys.stderr, self.tostring()
		
			return (-1, -1)

		i = rlen - 1
		while i >= 0:
			if up_seq[i].upper() == down_seq[i].upper():
				i -= 1
			else:
				break
		j = rlen
		while j < 2*rlen:
			if up_seq[j].upper() == down_seq[j].upper():
				j += 1
			else:
				break

		print up_seq.upper()[0:rlen], up_seq.lower()[rlen:]
		print down_seq.lower()[0:rlen], down_seq.upper()[rlen:]
		return (rlen-1-i, j-rlen)
				
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
	def check_bed_ann_list(self, bed_ann_list):
		flag = True
		if not bed_ann_list:
			flag = False
		else:
			if len(set([(a.gene_name, a.chr, a.strand) for a in bed_ann_list])) > 1:
				flag = False
			else:
				max_start = max([a.start for a in bed_ann_list])
				min_end = min([a.end for a in bed_ann_list])
				# Gene has two annotations more than 1Mb away from each other
				if max_start - min_end > 5000000:
					flag = False
			print >> sys.stderr, "Warnning: Input gene annotations include multiple chr, strand, or regions (5Mb away). Skipping current gene annotation."
			print >> sys.stderr, set([(a.gene_name, a.chr, a.strand) for a in bed_ann_list])

		return flag
	def load(self, bed_ann_list):
		if not self.check_bed_ann_list(bed_ann_list):
			self.gene_name = "NA"
			self.chr = "NA"
			self.start = -1
			self.end = -1
			self.strand = "NA"
			self.is_coding = False
			self.is_contradictory = True # contradictory gene annotation, will not be used
		else:
			self.gene_name = bed_ann_list[0].gene_name
			self.chr = bed_ann_list[0].chr
			self.strand= bed_ann_list[0].strand
			self.start = min([a.start for a in bed_ann_list])
			self.end = max([a.end for a in bed_ann_list])
			self.is_coding = True if "cds" in [a.type for a in bed_ann_list] else False
			self.is_contradictory = False
		'''
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
				print >> sys.stderr, "Warnning: GeneBed annotation includes multiple genes. Only the first annotation will be used."
				print >> sys.stderr, set([(a.gene_name, a.chr, a.strand) for a in bed_ann_list])
				#sys.exit(1)
			self.gene_name = bed_ann_list[0].gene_name
			self.chr = bed_ann_list[0].chr
			self.strand= bed_ann_list[0].strand

			self.start = min([a.start for a in filter(lambda x:x.chr==self.chr, bed_ann_list)])
			self.end = max([a.end for a in filter(lambda x:x.chr==self.chr, bed_ann_list)])
			#self.start = min([a.start for a in bed_ann_list])
			#self.end = max([a.end for a in bed_ann_list])
			#self.is_coding = True if "cds" in [a.type for a in bed_ann_list] else False
			self.is_coding = True if "cds" in [a.type for a in filter(lambda x:x.chr==self.chr, bed_ann_list)] else False
		'''		
	
	def overlap(self, interval2):
		if self.chr == interval2.chr and min(self.end, interval2.end) - max(self.start, interval2.start) > 0:
			return True
		else:
			return False
	def merge(self, interval2):
		if self.chr != interval2.chr or self.strand != interval2.strand:
			print >> sys.stderr, "Warning: intervals are on different chr/strand."
			print >> sys.stderr, self.gene_name, interval2.gene_name
			return self
		else:
			#new_interval = GeneIntervals("Merged_" + self.gene_name + "_" + interval2.gene_name, self.chr, min(self,start, interval2.start), max(self.end, interval2.end), self.strand, False)
			new_interval = copy.copy(interval2)
			new_interval.gene_name = "Merged_" + self.gene_name + "_" + interval2.gene_name
			new_interval.chr =  self.chr
			new_interval.start =  min(self.start, interval2.start)
			new_interval.end =  max(self.end, interval2.end)
			new_interval.strand = self.strand
			# for merged gene intervals is_coding has no sense
			new_interval.is_coding = False
			
			return new_interval			
# Deprecated, use GeneBed instead
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
	def print_coding_gene_idx(self):
		for key in self.__gene_name_idx_map:
			print key, self.__gene_name_idx_map[key]
	# whether a gene include coding exon (cds)
	def is_coding(self, gene_name):
		return self.__gene_intervals[gene_name].is_coding 
	# where a gene has contradictory annotations
	def is_contradictory(self, gene_name):
		return self.__gene_intervals[gene_name].is_contradictory
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
		# if pos is within 10bp window of any boundary, set close_to_boundary True
		t = 10
		matched_genes = []
		if not chr in self.__gene_starts:
			return matched_genes
		idx = bisect.bisect(self.__gene_starts[chr], pos)	
		while 0 < idx <= len(self.__gene_starts[chr]):
			#bpann is an GeneBed object
			bpann = self.__genes[chr][idx-1]
			#search within a limited region (default 1000000)
			if pos - bpann.start > self.__max_diff: 
				break

			if bpann.start <= pos <= bpann.end:
				# check if pos is on/close to current annotation's boundary, these are not gene annotations' attributes, need to reset every time
				bpann.is_on_boundary = False
				bpann.close_to_boundary = False
				if bpann.start == pos or bpann.end == pos:
					bpann.is_on_boundary = True
					bpann.close_to_boundary = True
				elif min(abs(bpann.start-pos), abs(bpann.end-pos)) < t:
					bpann.close_to_boundary = True

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
						if bpann.transcript_id == adjacent_bpann.transcript_id and adjacent_bpann.type == "cds" and abs(adjacent_bpann.idx - bpann.idx) <= 1:
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
						if bpann.transcript_id == adjacent_bpann.transcript_id and adjacent_bpann.type == "cds" and abs(adjacent_bpann.idx - bpann.idx) <= 1:
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
		
		
def cmp_fusion_breakpoints(bp1, bp2, diff):
	chr1, pos1, strand1 = bp1
	chr2, pos2, strand2 = bp2
			
	if chr1 != chr2 or strand1 != strand2:
		return False
	elif abs(pos1 - pos2) < diff:
		return True
	else:
		return False	
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
