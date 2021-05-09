import sys
import re

def del_emp(l):
	l2 = []
	for i in range(len(l)):
		if not l[i] == "":
			l2.append(l[i])
	return l2

gene_dict = {}
f = open(sys.argv[1])
for line in f:
	line = line.replace("\n","")
	line_l = line.split("\t")
#	print(line)
	if line_l[2] == "UNANNOTATED":
#		print(line)
		if not line_l[3] in gene_dict:
			gene_dict[line_l[3]] = []
			gene_dict[line_l[3]].append(line_l[5])
		else:
			gene_dict[line_l[3]].append(line_l[5])

#for k,v in gene_dict.items():
#	print(k,v,sep="\t")

exonnum_pos_dict = {}
f3 = open(sys.argv[3])
for line in f3:
	line = line.replace("\n","")
	line_l = line.split("\t")
#	print(line)
	exonnum_l = line_l[16].split(",")
	start_l = del_emp(line_l[9].split(","))
	end_l = del_emp(line_l[10].split(","))
	for i in range(len(exonnum_l)):
		if not line_l[12] +"/"+ exonnum_l[i] in exonnum_pos_dict:
			exonnum_pos_dict[line_l[12] +"/"+ exonnum_l[i]] = [start_l[i], end_l[i]]

#for k,v in exonnum_pos_dict.items():
#	print(k,v,sep="\t")

convert_dict = {}
f2 = open(sys.argv[2])
for line in f2:
	line = line.replace("\n","")
	line_l = line.split("\t")
#	print(line)
	if line.startswith("Variant_type"):
		continue
	info_l = line_l[0].split(";")
	gene_l = info_l[1].split("/")
	gene_name = ""
	for i in range(len(gene_l)):
		gene_l[i]
#		print(gene_l[i])
		if not gene_l[i] =="novel":
			gene_name = gene_l[i]
	if not gene_name in gene_dict:
#		print(">>>")
		continue
#	print(line)
#	print("gene_name: ", gene_name, sep="\t")
#	print("gene_dict[gene_name]: ", gene_dict[gene_name], sep="\t")	
#	print("info_l: ", info_l, sep="\t")	
	exonnum_l = del_emp(re.split('[*/,]',info_l[2]))
#	print("exonnum_l: ", exonnum_l, sep="\t")
	eva_len_l = info_l[3].split(",")
#	print("eva_len_l: " , eva_len_l, sep="\t")
	
	convert_exonnum_l = []
	exonnum_count = 0
	for i in range(len(eva_len_l)):
#		print(eva_len_l[i])
		if "chr" in eva_len_l[i]:
			convert_exonnum_l.append("novel")
		else:
			convert_exonnum_l.append(exonnum_l[exonnum_count])
			exonnum_count += 1
#	print("convert_exonnum_l: ", ",".join(convert_exonnum_l), sep="\t")
	convert_exonpos_l = []
	for i in range(len(convert_exonnum_l)):
		if not convert_exonnum_l[i] == "novel":
#			print(convert_exonnum_l[i], exonnum_pos_dict[gene_name+"/"+convert_exonnum_l[i]], sep="\t")
			convert_exonpos_l.append(exonnum_pos_dict[gene_name+"/"+convert_exonnum_l[i]][0])
			convert_exonpos_l.append(exonnum_pos_dict[gene_name+"/"+convert_exonnum_l[i]][1])
		else:
#			print(convert_exonnum_l[i])	
			convert_exonpos_l.append("novel")
			convert_exonpos_l.append("novel")
#	print("convert_exonpos_l: ", convert_exonpos_l, sep="\t")
	convert_exonpos_l2 = convert_exonpos_l
	del convert_exonpos_l2[0]
	del convert_exonpos_l2[-1] 
#	print("convert_exonpos_l2: ", ",".join(convert_exonpos_l2), sep="\t")
	
	
	for m_pos in gene_dict[gene_name]:
#		print(m_pos)#2,4,5,novel/s,k,k,chr4_73484454_73485677/*,*/*
		m_pos_l = m_pos.split("/")
#		print("m_pos_l: ", m_pos_l, sep="\t")
		m_pos_convert_exon_pos_l = []
		m_pos_exon_num_l = m_pos_l[0].split(",")
		for i in range(len(m_pos_exon_num_l)):
#			print(exonnum_pos_dict[gene_name+"/"+m_pos_exon_num_l[i]])
			if not m_pos_exon_num_l[i] == "novel":
				m_pos_convert_exon_pos_l.append(exonnum_pos_dict[gene_name+"/"+m_pos_exon_num_l[i]][0])
				m_pos_convert_exon_pos_l.append(exonnum_pos_dict[gene_name+"/"+m_pos_exon_num_l[i]][1])
			else:
				m_pos_convert_exon_pos_l.append("novel")
				m_pos_convert_exon_pos_l.append("novel")
#		print("m_pos_convert_exon_pos_l: ", ",".join(m_pos_convert_exon_pos_l), sep="\t")
		if ",".join(convert_exonpos_l2) in ",".join(m_pos_convert_exon_pos_l):
#			print("OVERRAP")
			if not m_pos in convert_dict:
				convert_dict[m_pos] = [",".join(convert_exonnum_l), info_l[3]]
#				print("convert_dict: ", convert_dict, sep="\t")
		
#	print()
	
#for k,v in convert_dict.items():
#	print(k,v,sep="\t")

f4 = open(sys.argv[4])
for line in f4:
	line = line.replace("\n","")
	line_l = line.split("\t")
#	print(line)
	if line_l[3] == "UNANNOTATED":
#		print(line)
		if line_l[6] in convert_dict:
#			print(line_l[0])
#			print(line_l[6])
#			print(convert_dict[line_l[6]])
			converted_l1_l = convert_dict[line_l[6]][1].split(",")
			
			converted_l2 = []
			for i in range(len(converted_l1_l)):
#				print(converted_l[i])
				if "chr" in converted_l1_l[i]:
					converted_l2.append(converted_l1_l[i])
				else:
					if converted_l1_l[i] == "known":
						converted_l2.append("k")
					elif converted_l1_l[i] == "short":
						converted_l2.append("s")
					elif converted_l1_l[i] == "long":
						converted_l2.append("l")
					elif converted_l1_l[i] == "slong":
						converted_l2.append("sl")
#			print("converted_l2: ", converted_l2, sep="\t")
			print(line_l[0], line_l[6], convert_dict[line_l[6]][0] +"/"+ ",".join(converted_l2) +"/*,*/*", sep="\t")		
		else:
			print(line_l[0], line_l[6], line_l[6], sep="\t")
#		print()
#	print()
	
