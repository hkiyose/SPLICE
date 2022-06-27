import sys
import re

def add_exonnum(f):
	def remove_duplicates(x):
		y=[]
#		dup_count = 0
		for i in x:
			if i not in y:
				y.append(i)
		return y

	def get_uniq_sort_dict(dict):
		for k,v in dict.items():
#			print(k,v,sep="\t")
#			print(len(v))
			uniq_v = remove_duplicates(v)
#			print("uniq_v: ", uniq_v, sep="\t")
			uniq_v_sort = sorted(uniq_v)
#			print("uniq_v_sort: ", uniq_v_sort, sep="\t")
			exon_dict_uniq[k] = uniq_v_sort
#			print("exon_dict_uniq[k]: ", len(exon_dict_uniq[k]), exon_dict_uniq[k], sep="\t")
		return exon_dict_uniq

	def get_exon_num(dict):
#		print(dict)
		for k,v in dict.items():
			exon_num = {}
#			print(k,v,sep="\t")
			for i in range(len(v)):
#				print(k, v[i][0], v[i][1])
				exon_pos = v[i][0],v[i][1]
				exon_num[exon_pos] = i+1
				exon_dict_uniq_num[k] = exon_num
		return exon_dict_uniq_num

	f1 = open(f)
	f2 = open(f)
	gene_dict = {}
	exon_dict = {}
	exon_dict_uniq = {}
	exon_dict_uniq_num = {}
	pre_gene_name = ""
	gene_list = []
	uniq_dict = {}

	for line in f1:
		line = line.replace("\n","")
		line_l = line.split("\t")
		if line.startswith("#"):
			continue
		exon_start = line_l[9].split(",")
		exon_end = line_l[10].split(",")
		del exon_start[-1]
		del exon_end[-1]
	
		for i in range(len(exon_start)):
			exon_dict.setdefault(line_l[12], []).append([int(exon_start[i]),int(exon_end[i])])
	#print(exon_dict)
	exon_uniq_sort_dict = get_uniq_sort_dict(exon_dict)
	#print(exon_uniq_sort_dict)
	exon_uniq_sort_num_dict = get_exon_num(exon_uniq_sort_dict)
	#print(exon_uniq_sort_num_dict)

	for line in f2:
		line = line.replace("\n","")
		line_l = line.split("\t")
		if line.startswith("#"):
			continue
		exon_start = line_l[9].split(",")
		exon_end = line_l[10].split(",")
		del exon_start[-1]
		del exon_end[-1]
	
		exon_num_list = []
		for i in range(len(exon_start)):
			exon_pos = int(exon_start[i]),int(exon_end[i])
			exon_num_list.append(exon_pos)
#	    	print(exon_num_list)
		exon_num_iso_list = []
		exon_num_iso_str = ""
		for exon in exon_num_list:
			exon_num_iso_list.append(int(exon_uniq_sort_num_dict[str(line_l[12])][exon]))
#	    print(exon_num_iso_list)
		mapped_list = map(str, exon_num_iso_list)
		mapped_list_str = ",".join(mapped_list)
#	    print(mapped_list_str)
		print(line+"\t"+mapped_list_str)


def rm_conjoined(f):
	f1 = open(f)
	f2 = open(f)
	gene_l = []
	for line in f1:
		line = line.replace("\n","")
		line_l = line.split("\t")
#		print(line)
		if not line_l[12] in gene_l and not "-" in line_l[12]:
			gene_l.append(line_l[12])
#			print(gene_l)
	for line in f2:
		line = line.replace("\n","")
		line_l = line.split("\t")
		eva_cg = False	
		if "-" in line_l[12]:
#			print(">>> - in gene_name", line_l[12], line, sep="\t")
			gene_l2 = line_l[12].split("-")
			if gene_l2[0] == "AS":
				print(line)
				continue
			else:	
				for i in range(len(gene_l)):
					if line_l[12] == gene_l2[0] + "-" + gene_l[i]:
#						print(">>> conjoined_gene")
#						print(line)
						eva_cg = True	
		if eva_cg == False:
			print(line)

def ref_convert(f):
	f1 = open(f)
	for line in f1:
		line = line.replace("\n","")
		line_l = line.split("\t")
#		print(line_l[12])
		if "/" in line_l[12]:
			gene_l = line_l[12].split("/")
			line_l[12] = gene_l[0]+gene_l[1]
			print("\t".join(line_l))
		else:
			print("\t".join(line_l))


def convert_fa(f):
	seq = ""
	f1 = open(f)
	for line in f1:
		line = line.replace("\n","")
		line_l = line.split("\t")
#		print(line)
		if line.startswith(">"):
			if not seq == "":
				print(seq)
				seq = ""
			print(line)
		else:
			seq += line
#			print(seq)

def make_ref_seq_cp(f, f_2):
	def del_terminal_index(l):
		if l[-1] == "":
			del l[-1]
		return l

	f1 = open(f)
	chr_seq_dict = {}
	for line in f1:
		line = line.replace("\n", "")
		line_l = line.split()
		if line.startswith(">"):
			if line.startswith(">NC_"):
				chr2 = line_l[0][1:]
			else:
				chr2 = line[1:]
			chr_seq_dict[chr2] = ""
		else:
			chr_seq_dict[chr2] = line
	
			f2 = open(f_2)
			seq = ""
			for line in f2:
				line = line.replace("\n", "")
				line_l = line.split("\t")
				start_pos_l = del_terminal_index(line_l[9].split(","))
				end_pos_l = del_terminal_index(line_l[10].split(","))
				if line_l[2] in chr_seq_dict:
					for i in range(len(start_pos_l)):
						seq += chr_seq_dict[line_l[2]][int(start_pos_l[i]):int(end_pos_l[i])]
					if not seq == "":
						print(">" + line_l[1] + "/" + line_l[12])
						print(seq)
				seq = ""
			chr_seq_dict = {}

