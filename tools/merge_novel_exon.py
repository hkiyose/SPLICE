import sys
import os
import re
from glob import glob
import math

def merge_novel_exon(d):

	dir_name = str(d)
	file_l = []
	for file in glob(dir_name + '/*.novel_exon'):
	#	print(file)
		file_l.append(file)

	variant_read_count_dict = {}#variant_read_count_dict[k] = [cancer_read_count, liver_read_count]
	variant_error_rate_dict = {}
	read_count_l = [0]*len(file_l)
	sample_name_l = []
	file_count = -1
	for file in sorted(file_l):
		file_count += 1
		f = open(file)
		basename = os.path.basename(file)
	#	print("basename: ", basename, sep="\t")
		basename_l = basename.split("_")
	#	print(basename_l[0], basename_l[1], sep="\t")	
		
		sample_name = basename_l[0]
	#	print("sample_name: ", sample_name, sep="\t")
		sample_name_l.append(sample_name)
		
		for line in f:
			line = line.replace("\n", "")
			line_l = line.split("\t")	
			read_count_l[file_count] += int(line_l[0])
			if line_l[0] == "1":
				continue
	#		print(line)	

			k = line_l[1] +";"+ line_l[2] +";"+ line_l[3] +";"+ line_l[4]
	#		print(k)
			
			if not k in variant_read_count_dict:
				variant_read_count_dict[k] = ["0"]*len(file_l)
				
	#			print(variant_read_count_dict[k])
			
			variant_read_count_dict[k][file_count] = line_l[0]
	#		print("sample_name_l: ", sample_name_l, sep="\t")
	#		print("variant_read_count_dict[k]: ", variant_read_count_dict[k], sep="\t")		
	#		print()

	print("Variant_type","\t".join(sample_name_l),sep="\t")
	for k,v in variant_read_count_dict.items():
		print(k,"\t".join(v),sep="\t")


def merge_novel_exon2(f):

	novel_exon_dict = {}
	f1 = open(f)
	for line in f1:
		line = line.replace("\n","")
		line_l = line.split("\t")
	#	print(line)
		if line.startswith("Variant"):
			continue
		info_l = line_l[0].split(";")
	#	print("info_l: ", info_l, sep="\t")
		chr_pos_l = info_l[3].split(",")
	#	print("chr_pos_l: ", chr_pos_l, sep="\t")

		for i in range(len(chr_pos_l)):
			if not chr_pos_l[i] == "known" and not chr_pos_l[i] == "short" and not chr_pos_l[i] == "long" and not chr_pos_l[i] == "slong":
				chr_pos_l2 = chr_pos_l[i].split("_")
				if len(chr_pos_l2) == 5:
	#				print("5!: ", line, sep="\t")
					k = chr_pos_l2[0]+"_"+chr_pos_l2[1]+"_"+chr_pos_l2[2], chr_pos_l2[3], chr_pos_l2[4]
					if k in novel_exon_dict:
						novel_exon_dict[k] += 1
					else:	
						novel_exon_dict[k] = 1
				elif len(chr_pos_l2) == 4:
	#				print("4!: ", line, sep="\t")
					k = chr_pos_l2[0]+"_"+chr_pos_l2[1], chr_pos_l2[2], chr_pos_l2[3]
					if k in novel_exon_dict:
						novel_exon_dict[k] += 1
					else:
						novel_exon_dict[k] = 1
				else:#elif len(chr_pos_l2) == 3:
	#				print("3!!!: ", line, sep="\t")
					k = chr_pos_l2[0], chr_pos_l2[1], chr_pos_l2[2]
					if k in novel_exon_dict:
						novel_exon_dict[k] += 1
					else:
						novel_exon_dict[k] = 1			
	#	print()

	for k,v in sorted(novel_exon_dict.items(), key=lambda x: (x[0][0], int(x[0][-2]), int(x[0][-1]))):
		print(k[0], k[1], k[2], v, sep="\t")


def merge_novel_exon3(f):

	merged_pos = []
	f1 = open(f)
	num_lines = sum(1 for line in open(f))
	#print(num_lines)
	line_count = 0
	for line in f1:
		line_count += 1
		line = line.replace("\n","")
		line_l = line.split("\t")
	#	print("line: ", line, sep="\t")
		
		if len(merged_pos) == 0:
			merged_pos = [line_l[0], line_l[1], line_l[2]]
	#		print("merged_pos: " ,merged_pos ,sep="\t")
	#		print()
			continue
		
		if int(line_l[1]) > int(merged_pos[2]) or line_l[0] != merged_pos[0]:
	#		print("!")
			print("\t".join(merged_pos))
			merged_pos = [line_l[0], line_l[1], line_l[2]]
	#		print("merged_pos: " , merged_pos, sep="\t")
			if num_lines == line_count:
				print("\t".join(merged_pos))
		else:
			novel_pos_in_merged_set = set([x for x in range(int(merged_pos[1]), int(merged_pos[2])+1)])
			novel_pos_in_line_set = set([x for x in range(int(line_l[1]), int(line_l[2])+1)])
			and_set = novel_pos_in_merged_set & novel_pos_in_line_set
	#		print("novel_pos_in_merged_set: ", len(novel_pos_in_merged_set), novel_pos_in_merged_set, sep="\t")
	#		print("novel_pos_in_line_set: ", len(novel_pos_in_line_set), novel_pos_in_line_set, sep="\t")
	#		print("and_set_count: ", len(and_set), and_set, sep="\t")
			if len(and_set) > 0:	
	#			print("overlap")
				if int(line_l[2]) > int(merged_pos[2]):
					merged_pos[2] = line_l[2]
	#				print("merged_pos: " , merged_pos, sep="\t")
			if num_lines == line_count:
				print("\t".join(merged_pos))
	#	print()	
	

def merge_novel_exon4(f, f_2):

	chr_pos_d = {}
	f1 = open(f)
	for line in f1:
		line = line.replace("\n","")
		line_l = line.split("\t")
	#	print(line)
		v = [line_l[1], line_l[2]]
		if line_l[0] in chr_pos_d:
			chr_pos_d[line_l[0]].append(v)	
		else:
			chr_pos_d[line_l[0]] = []
			chr_pos_d[line_l[0]].append(v)

	#for k,v in chr_pos_d.items():
	#	print(k,v,sep="\t")

	novel_pos_dict = {}
	f2 = open(f_2)
	for line in f2:
		line = line.replace("\n","")
		line_l = line.split("\t")
	#	print(line)
		if line.startswith("Vari"):
			print(line)
			continue
		
		info_l = line_l[0].split(";")
	#	print("info_l: ", info_l, sep="\t")
		
		novel_pos_l = info_l[-1].split(",")
	#	print("novel_pos_l: ", novel_pos_l, sep="\t")

		converted_novel_pos_l = []
		
		for i in range(len(novel_pos_l)):
	#		print(novel_pos_l[i])
			if not "_" in novel_pos_l[i]:
				converted_novel_pos_l.append(novel_pos_l[i])
			else:
	#			print(novel_pos_l[i])
				pos_l = novel_pos_l[i].split("_")
	#			print("pos_l: ", pos_l, sep="\t")	
				chr_s = ""
				if len(pos_l) == 5:	
					chr_s = "_".join(pos_l[0:3])
				elif len(pos_l) == 4:
					chr_s = "_".join(pos_l[0:2])
				else:
					chr_s = "_".join(pos_l[0:1])	
	#			print("chr_s: ", chr_s,sep="\t")
		
				for j in range(len(chr_pos_d[chr_s])):
	#				print(chr_pos_d[chr_s][j][0], chr_pos_d[chr_s][j][1], sep="\t")
					if int(chr_pos_d[chr_s][j][0]) <= int(pos_l[-2]) <= int(chr_pos_d[chr_s][j][1]):
	#					print(chr_pos_d[chr_s][j][0], chr_pos_d[chr_s][j][1], sep="\t")
						converted_novel_pos_l.append(chr_s +"_"+ chr_pos_d[chr_s][j][0] +"_"+ chr_pos_d[chr_s][j][1])
	#	print("converted_novel_pos_l:" ,converted_novel_pos_l, sep="\t")
		
		converted_info = ";".join(info_l[0:3]) +";"+ ",".join(converted_novel_pos_l)
	#	print("converted_info: ", converted_info, sep="\t")
		
		if converted_info in novel_pos_dict:
			for i in range(len(line_l[1:])):
	#			print(novel_pos_dict[converted_info][i], line_l[i+1])
				novel_pos_dict[converted_info][i] = str(int(novel_pos_dict[converted_info][i]) + int(line_l[i+1]))
	#			print(novel_pos_dict[converted_info][i], line_l[i+1])
	#		print("novel_pos_dict[converted_info]: ", novel_pos_dict[converted_info], sep="\t")
		else:
			novel_pos_dict[converted_info] = line_l[1:]
	#		print("novel_pos_dict[converted_info]: ", novel_pos_dict[converted_info], sep="\t")
	#	print()

	for k,v in novel_pos_dict.items():
		print(k,"\t".join(v), sep="\t")

def merge_novel_exon5(f):

	f1 = open(f)
	for line in f1:
		line = line.replace("\n","")
		line_l = line.split("\t")
		if line.startswith("Variant"):
			continue
	#	print(line)
		info_l = line_l[0].split(";")
	#	print("info_l: ", info_l, sep="\t")
		
		annot_gene_l = info_l[1].split("/")
		annot_gene = "" 
		for i in range(len(annot_gene_l)):
	#		print(annot_gene_l[i])
			if annot_gene_l[i] == "novel":
				continue
			else:
				if annot_gene == "":
					annot_gene = annot_gene_l[i]
				else:
					annot_gene += "/"+annot_gene_l[i]
		
		if annot_gene == "":
			annot_gene = "NOVEL"
	#	print("annot_gene: ", annot_gene, sep="\t")
		
		eva_exon_len_l = info_l[3].split(",")
		convert_eva_exon_len_l = []
		for i in range(len(eva_exon_len_l)):
			if eva_exon_len_l[i] == "known":
				convert_eva_exon_len_l.append("k")
			elif eva_exon_len_l[i] == "short":
				convert_eva_exon_len_l.append("s")
			elif eva_exon_len_l[i] == "long":
				convert_eva_exon_len_l.append("l")
			elif eva_exon_len_l[i] == "slong":
				convert_eva_exon_len_l.append("sl")
			else:
				convert_eva_exon_len_l.append(eva_exon_len_l[i])

		exon_num_l = re.split('[,/]', info_l[2])
		convert_exon_num_l = []
	#	print("exon_num_l: ", exon_num_l, sep="\t")
		exon_num_l2 = []
		for i in range(len(exon_num_l)):
			if exon_num_l[i] == "*":
				continue
			else:
				exon_num_l2.append(exon_num_l[i])
	#	print("exon_num_l2: ", exon_num_l2, sep="\t")
		
		count = 0
		for i in range(len(convert_eva_exon_len_l)):
			if "_" in convert_eva_exon_len_l[i]:
				convert_exon_num_l.append("novel")		
			else:
				convert_exon_num_l.append(exon_num_l2[count])
				count += 1
				
	#	print("convert_exon_num_l: ", convert_exon_num_l, sep="\t")
		
		print(annot_gene, "NOVEL", ",".join(convert_exon_num_l) +"/"+ ",".join(convert_eva_exon_len_l) +"/*,*/*", "\t".join(line_l[1:]), sep="\t")	
		
	#	print()


def merge_novel_exon6(f):

	pre_gene = ""
	f1 = open(f)
	for line in f1:
		line = line.replace("\n","")
		line_l = line.split("\t")
	#	print(line)
		if line_l[3] != pre_gene:
			pre_gene = line_l[3]
			gene_variant_l =[]
			variant_dict = {}
			gene_variant_l.append(line_l[5])
			start_pos_l = line_l[7].split(",")
			end_pos_l = line_l[8].split(",")
	#		print("start_pos_l: ", start_pos_l, sep="\t")	
	#		print("end_pos_l: ", end_pos_l, sep="\t")
			full_pos_l = []
			for i in range(len(start_pos_l)):
	#			print(start_pos_l[i], end_pos_l[i], sep="\t")
				full_pos_l.append(start_pos_l[i])
				full_pos_l.append(end_pos_l[i])
	#		print("full_pos_l: ", len(full_pos_l), full_pos_l, sep="\t")
			full_pos_s = ",".join(full_pos_l)
	#		print("full_pos_s: ", full_pos_s, sep="\t")
			variant_dict[line_l[5]] = full_pos_s	
	#		print("variant_dict[line_l[5]]: ", variant_dict[line_l[5]], sep="\t")
			print("\t".join(line_l[0:6]), "\t".join(line_l[9:]), sep="\t")
	#		print()	
			continue
			
		else:
			start_pos_l = line_l[7].split(",")	
			end_pos_l = line_l[8].split(",")
	#		print("start_pos_l: ", start_pos_l, sep="\t")
	#		print("end_pos_l: ", end_pos_l, sep="\t")
			full_pos_l = []
			for i in range(len(start_pos_l)):
	#			print(start_pos_l[i], end_pos_l[i], sep="\t")
				full_pos_l.append(start_pos_l[i])
				full_pos_l.append(end_pos_l[i])
	#		print("full_pos_l: ", len(full_pos_l), full_pos_l, sep="\t")	
			full_pos_s = ",".join(full_pos_l)
	#		print("full_pos_s: ", full_pos_s, sep="\t")	
			
			del full_pos_l[0]
			del full_pos_l[-1]
			full_pos_s2 = ",".join(full_pos_l)
	#		print("full_pos_s2: ", full_pos_s2, sep="\t") 
				
			eva_overlap = False	
			for i in range(len(gene_variant_l)):
	#			print(">>> gene_variant_l[i]: ", gene_variant_l[i], sep="\t")
				pre_gene_s = variant_dict[gene_variant_l[i]]
	#			print("pre_gene_s: ", pre_gene_s, sep="\t")
				
				if eva_overlap == True:
					continue
		
				if full_pos_s2 in pre_gene_s:
	#				print("OVERLAP", "MERGE")
					print("\t".join(line_l[0:4]), "PARTIAL", gene_variant_l[i], "\t".join(line_l[9:]), sep="\t")
					eva_overlap = True
	#				print("--------------------")	
	#			else:
	#				print("NOT OVERLAP")		
	#				print("--------------------")	
		
			if eva_overlap == False:
				gene_variant_l.append(line_l[5])
				variant_dict[line_l[5]] = full_pos_s
	#			print("gene_variant_l: ", len(gene_variant_l), gene_variant_l, sep="\t")
	#			print("variant_dict[line_l[5]]: ", variant_dict[line_l[5]], sep="\t")
				print("\t".join(line_l[0:6]), "\t".join(line_l[9:]), sep="\t")	
	#		print()	 
			
	#	print()
		
