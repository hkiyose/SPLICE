import sys
import os
import re
from glob import glob
import math

def merge_table(d):

	dir_name = str(d)
	file_l = []
	for file in glob(dir_name + '/*.annot'):
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
	#		print(line)	
			read_count_l[file_count] += int(line_l[0])
			k = line_l[2] +"/"+ line_l[3] +"/"+ line_l[4] +"/"+ line_l[5] +"/"+ line_l[6] +"/"+ line_l[7] +"/"+ line_l[10]	
			if not k in variant_read_count_dict:
				variant_read_count_dict[k] = ["0"]*len(file_l)
				variant_error_rate_dict[k] = ["0"]*len(file_l)	
	#			print(variant_read_count_dict[k])
	#			print(variant_read_rate_dict[k])
			
			variant_read_count_dict[k][file_count] = line_l[0]
			variant_error_rate_dict[k][file_count] = line_l[8]
	#		print("sample_name_l: ", sample_name_l, sep="\t")
	#		print("variant_read_count_dict: ", variant_read_count_dict, sep="\t")		
	#		print("variant_read_rate_dict: ", variant_read_rate_dict, sep="\t")
	#		print()
	#	print()
				
	sample_name_read_rate_l = []
	sample_name_error_rate_l = []
	sample_name_fc_l = []
	sample_name_log2_fc_l = []
	for i in range(len(sample_name_l)):
		sample_name_read_rate_l.append(sample_name_l[i] + "_read/1Mread")
		sample_name_error_rate_l.append(sample_name_l[i] + "_error_rate")
		
		sample_name_l_l = sample_name_l[i].split("_")
		if not sample_name_l_l[0] in sample_name_fc_l:
			if not sample_name_l_l[0] + "_FC(Cancer/Liver)" in sample_name_fc_l:
				sample_name_fc_l.append(sample_name_l_l[0] + "_FC(Cancer/Liver)")
				sample_name_log2_fc_l.append(sample_name_l_l[0] + "_log2(FC(Cancer/Liver))")

	print("Log2_FC_ave" , "variant_info", "\t".join(sample_name_l), "\t".join(sample_name_read_rate_l), "\t".join(sample_name_fc_l), "\t".join(sample_name_log2_fc_l), "\t".join(sample_name_error_rate_l), sep="\t")

	#print(read_count_l)

	for k,v in variant_read_count_dict.items():
	#	print(len(read_count_l), read_count_l, sep="\t")
	#	print(len(v), v, sep="\t")
		read_rate_l = [""]*len(file_l)
		for i in range(len(read_count_l)):
			if not v[i] == "0":
				read_rate_l[i] = str(int(v[i])/int(read_count_l[i])*1000000)
			else:
				read_rate_l[i] = "0"
	#	print(len(read_rate_l), read_rate_l, sep="\t")

	#	print(str(log2_fc_ave), k, "\t".join(v), "\t".join(read_rate_l), "\t".join(fc_l), "\t".join(log2_fc_l), "\t".join(variant_error_rate_dict[k]), sep="\t")
		print("-", k, "\t".join(v), "\t".join(read_rate_l), "-", "-", "\t".join(variant_error_rate_dict[k]), sep="\t")
	#	print()


def merge_known(f):

	transcript_read_count_dict = {}
	f1 = open(f)
	for line in f1:
		line = line.replace("\n","")
		line_l = line.split("\t")
		if line_l[1] == "NOVEL":
			continue
#		print(line)
		k = "/".join(line_l[0:6])
#		print("k: ", k, sep="\t")
		
		if not k in transcript_read_count_dict:
			transcript_read_count_dict[k] = []
			for reads in line_l[6:]:
				transcript_read_count_dict[k].append(int(float(reads)))
#			print("transcript_read_count_dict[k]: ", transcript_read_count_dict[k], sep="\t")
		else:
			for i in range(len(line_l[6:])):
#				print(line_l[i+6])
				transcript_read_count_dict[k][i] += int(float(line_l[i+6]))  
#			print("transcript_read_count_dict[k]: ", transcript_read_count_dict[k], sep="\t")	
#		print()
	
	for k,v in transcript_read_count_dict.items():
		v_s = []
		for i in range(len(v)):
			v_s.append(str(v[i]))
		k_l = k.split("/")
		print("\t".join(k_l), "\t".join(v_s), sep="\t")


def merge_len_comb(f, f_2, str3):

	def del_term_index(l):
		if l[-1] == "":
			del l[-1]
		return l

	def del_empty_index(l):
		l2 = []
		for i in range(len(l)):
			if not l[i] == "" and not l[i] == "*":
				l2.append(l[i])
		return l2

	gene_pos_dict = {}#gene_pos_dict[gene/exon_num] = [start_pos, end_pos]
	f1 = open(f)
	for line in f1:
		line = line.replace("\n", "")
		line_l = line.split("\t")
	#	print(line)
		start_pos_l = del_term_index(line_l[9].split(","))
		end_pos_l = del_term_index(line_l[10].split(","))
		exon_num_l = line_l[-1].split(",")
	#	print("start_pos_l: ", start_pos_l, sep="\t")
	#	print("end_pos_l: ", end_pos_l, sep="\t")
	#	print("exon_num_l: ", exon_num_l, sep="\t")
		for i in range(len(start_pos_l)):
			k = line_l[12] +"/"+ exon_num_l[i]
			if not k in gene_pos_dict:
				gene_pos_dict[k] = [start_pos_l[i], end_pos_l[i]]

	f2 = open(f_2)
	for line in f2:
		line = line.replace("\n","")
		line_l = line.split("\t")
		if not line_l[1] == "NOVEL":
			continue
		if not line_l[2] == str(str3):
			continue
	#	print(line)
		info_l = line_l[5].split("/")
		exon_num_l = del_empty_index(info_l[0].split(","))
		merged_exon_pos = info_l[2].split(",")
		gap_l = del_empty_index(info_l[3].split(","))
	#	print("info_l: ", info_l, sep="\t")
	#	print("exon_num_l: ", exon_num_l, sep="\t")
	#	print("merged_exon_pos: ", merged_exon_pos, sep="\t")
	#	print("gap_l: ", gap_l, sep="\t")
		
		pos_l = []	
		pos_l.append(merged_exon_pos[0])
		for i in range(len(exon_num_l)):
	#		print(gene_pos_dict[line_l[3]+"/"+exon_num_l[i]])
			pos_l.append(gene_pos_dict[line_l[3]+"/"+exon_num_l[i]][0])
			pos_l.append(gene_pos_dict[line_l[3]+"/"+exon_num_l[i]][1])
		pos_l.append(merged_exon_pos[1])
		
		pos_l2 = []
		start_pos_l = []
		end_pos_l = []	
	#	print("pos_l: ", pos_l, sep="\t")
		for i in range(len(pos_l)):
	#		print("pos: ", pos_l[i], "gap: ", gap_l[i], sep="\t")
			if i == 0 or i == len(pos_l)-1:
	#			print("nochange")
				pos_l2.append(pos_l[i])
			else:
				pos_l2.append(str(int(pos_l[i])+int(gap_l[i])))
				if i%2 == 1:
					start_pos_l.append(str(int(pos_l[i])+int(gap_l[i])))
				else:
					end_pos_l.append(str(int(pos_l[i])+int(gap_l[i])))
	#	print("pos_l2: ", pos_l2, sep="\t")
	#	print("start_pos_l: ", start_pos_l, sep="\t")
	#	print("end_pos_l: ", end_pos_l, sep="\t")
		
		transcript_length = 0
		for i in range(len(start_pos_l)):
			transcript_length += int(end_pos_l[i]) - int(start_pos_l[i]) + 1
	#	print("transcript_length: ", transcript_length, sep="\t")	
		print("\t".join(line_l[0:6]), transcript_length, ",".join(pos_l2), "\t".join(line_l[6:]), sep="\t")
	#	print()

	
def merge_len_comb2(f):

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
			pos_s = line_l[7]
	#		print("pos_s: ", pos_s, sep="\t")
			variant_dict[line_l[5]] = pos_s
	#		print("variant_dict[line_l[5]]: ", variant_dict[line_l[5]], sep="\t")
			print("\t".join(line_l[0:6]), "\t".join(line_l[8:]), sep="\t")
	#		print()
			continue
		else:
			pos_s = line_l[7]
	#		print("pos_s: ", pos_s, sep="\t")
			eva_overlap = False
			for i in range(len(gene_variant_l)):
	#			print(">>> gene_variant_l[i]: ", gene_variant_l[i], sep="\t")
				pre_gene_s = variant_dict[gene_variant_l[i]]
	#			print("pre_gene_s: ", pre_gene_s, sep="\t")
				
				if eva_overlap == True:
					continue
				
				if pos_s in pre_gene_s:
	#				print("OVERLAP", "MERGE")
					print("\t".join(line_l[0:4]), "PARTIAL", gene_variant_l[i], "\t".join(line_l[8:]), sep="\t")	
					eva_overlap = True
	#				print("--------------------")
	#			else:
	#				print("NOT OVERLAP")
	#				print("--------------------")
			if eva_overlap == False:
				gene_variant_l.append(line_l[5])
				variant_dict[line_l[5]] = pos_s
	#			print("gene_variant_l: ", len(gene_variant_l), gene_variant_l, sep="\t")
	#			print("variant_dict[line_l[5]]: ", variant_dict[line_l[5]], sep="\t")
				print("\t".join(line_l[0:6]), "\t".join(line_l[8:]), sep="\t")
	#		print()
	#	print()


def merge_len_comb3(f):
	
	transcript_read_count_dict = {}
	f1 = open(f)
	for line in f1:
		line = line.replace("\n","")
		line_l = line.split("\t")
	#	print(line)
		k = "/".join(line_l[0:6])
	#	print("k: ", k, sep="\t")

		if not k in transcript_read_count_dict:
			transcript_read_count_dict[k] = []
			for reads in line_l[6:]:
				transcript_read_count_dict[k].append(int(float(reads)))
	#		print("transcript_read_count_dict[k]: ", transcript_read_count_dict[k], sep="\t")
		else:
			for i in range(len(line_l[6:])):
				transcript_read_count_dict[k][i] += int(float(line_l[i+6]))
	#		print("MERGE")
	#		print("transcript_read_count_dict[k]: ", transcript_read_count_dict[k], sep="\t")
	#	print()

	for k,v in transcript_read_count_dict.items():
		v_s = []
		for i in range(len(v)):
			v_s.append(str(v[i]))
		k_l = k.split("/")
		print("\t".join(k_l[0:5]), "/".join(k_l[5:]), "\t".join(v_s), sep="\t")	


def merge_3utr(f):

	UTR3_merge_dict = {}
	f1 = open(f)
	for line in f1:
		line = line.replace("\n","")
		line_l = line.split("\t")
#		print(line)
		if line.startswith("CodingType"):
			print(line)
			continue
		transcript_info = ""
		if "3UTR" in line_l[5]:
			transcript_info = ";".join(line_l[0:5])+";"+"-"
		else:
			transcript_info = ";".join(line_l[0:6])
#		print("transcript_info: " ,transcript_info, sep="\t")
		
		if not transcript_info in UTR3_merge_dict:
			UTR3_merge_dict[transcript_info] = []
			for reads in line_l[10:]:
				UTR3_merge_dict[transcript_info].append(int(reads))
		else:
#			print("MERGE!")
			for i in range(len(line_l[10:])):
				UTR3_merge_dict[transcript_info][i] += int(line_l[i+10])
#		print(UTR3_merge_dict[transcript_info])
#		print()
	
	for k,v in UTR3_merge_dict.items():
		k_l = k.split(";")
		v_s = []
		for i in range(len(v)):
			v_s.append(str(v[i]))
		print("\t".join(k_l), "*", "*", "*", "*", "\t".join(v_s), sep="\t")
	
