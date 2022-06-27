import sys
import re

def novel_exon_filter(f, f_2, f_3):

	protein_coding_dict = {}
	f1 = open(f)
	for line in f1:
		line = line.replace("\n","")
		line_l = line.split("\t")
	#	print(line)
		if line_l[1] == "protein_coding":
			protein_coding_dict[line_l[0]] = 0

	f2 = open(f_2)
	for line in f2:
		line = line.replace("\n","")
		line_l = line.split("\t")
		if not line_l[0] in protein_coding_dict and line_l[1] == "protein_coding":
			protein_coding_dict[line_l[0]] = 0

	#for k,v in protein_coding_dict.items():
	#	print(k,v,sep="\t")

	f3 = open(f_3)
	for line in f3:
		line = line.replace("\n","")
		line_l = line.split("\t")
	#	print(line)
		if line_l[0] == "NOVEL":
			continue
	#	print(line)
		info_l = line_l[2].split("/")
		if line_l[0] in protein_coding_dict:
			print("CODING", "NOVEL", "UNANNOTATED", line_l[0], "-", line_l[2], "\t".join(line_l[3:]), sep="\t")
		else:
			print("NON-CODING", "NOVEL", "UNANNOTATED", line_l[0], "-", line_l[2], "\t".join(line_l[3:]), sep="\t")


def novel_exon_filter2(f, f_2):

	def del_term_index(l):
		if l[-1] == "":
			del l[-1]
		return l

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

	#for k,v in gene_pos_dict.items():
	#	print(k,v,sep="\t")	

	pre_gene = ""
	f2 = open(f_2)
	for line in f2:
		line = line.replace("\n", "")
		line_l = line.split("\t")
	#	print(line)
		if "/" in line_l[3]:	
	#		print()
			continue
		if "Log2_FC_ave" in line:
			continue
	#	print(line)
		info_l = line_l[5].split("/")
		exon_num_l = info_l[0].split(",")
	#	print("info_l: ", info_l, sep="\t")
		eva_exon_len_l = info_l[1].split(",")
	#	print("info_l: ", info_l, sep="\t")
	#	print("exon_num_l: ", exon_num_l, sep="\t")
	#	print("eva_exon_len_l: ", eva_exon_len_l, sep="\t")

		start_pos_l2 = []
		end_pos_l2 = []
		
		transcript_length = 0
		for i in range(len(exon_num_l)):
			if exon_num_l[i] == "novel":
				continue
	#		print(exon_num_l[i], gene_pos_dict[line_l[3]+"/"+exon_num_l[i]], sep="\t")
			transcript_length += int(gene_pos_dict[line_l[3]+"/"+exon_num_l[i]][1]) - int(gene_pos_dict[line_l[3]+"/"+exon_num_l[i]][0]) + 1
			start_pos_l2.append(int(gene_pos_dict[line_l[3]+"/"+exon_num_l[i]][0]))
			end_pos_l2.append(int(gene_pos_dict[line_l[3]+"/"+exon_num_l[i]][1]))
		
		for i in range(len(eva_exon_len_l)):
			if not "_" in eva_exon_len_l[i]:
				continue
	#		print(eva_exon_len_l[i])
			eva_exon_len_l_l = eva_exon_len_l[i].split("_")
	#		print(eva_exon_len_l_l[-2], eva_exon_len_l_l[-1], sep="\t")
			transcript_length += int(eva_exon_len_l_l[-1]) - int(eva_exon_len_l_l[-2]) + 1
			start_pos_l2.append(int(eva_exon_len_l_l[-2]))
			end_pos_l2.append(int(eva_exon_len_l_l[-1]))	
		
	#	print("start_pos_l2: ", sorted(start_pos_l2), sep="\t")	
	#	print("end_pos_l2: ", sorted(end_pos_l2), sep="\t")
	#	print("transcript_length: ", transcript_length, sep="\t")
		
		start_pos_l3 = []
		end_pos_l3 = []
		for pos in sorted(start_pos_l2):
			start_pos_l3.append(str(pos))	
		for pos in sorted(end_pos_l2):
			end_pos_l3.append(str(pos))
		
		print("\t".join(line_l[0:6]), transcript_length, ",".join(start_pos_l3), ",".join(end_pos_l3), "\t".join(line_l[6:]), sep="\t")
			
	#	print()


def novel_exon_sort(f):

	transcript_support_dict = {}
	f1 = open(f)
	for line in f1:
		line = line.replace("\n","")
		line_l = line.split("\t")
#		print(line)
		if line.startswith("Variant_type"):
			print(line)
			continue
		info_l = line_l[0].split(";")
		if info_l[1] == "novel":
			continue
		gene_l = info_l[1].split("/")
		if len(gene_l) > 2:
			continue
			
		total_reads = 0 
		for i in range(len(line_l)):
			if i == 0:
				continue
			total_reads += int(line_l[i])
#		print(line)
#		print(total_reads)
		transcript_support_dict[line] = total_reads
#		print()
		
	for k,v in sorted(transcript_support_dict.items(), key=lambda x:x[1], reverse=True):
		print(k)


def novel_exon_convert(f):

	splicing_variant_name_dict = {}
	f1 = open(f)
	for line in f1:
		line = line.replace("\n","")
		line_l = line.split("\t")
#		print(line)
		if line.startswith("Gene"):
			continue
		sample_name = "/".join(line_l[1:5])
#		print(sample_name)
		if sample_name in splicing_variant_name_dict:
			splicing_variant_name_dict[sample_name] += 1
		else:
			splicing_variant_name_dict[sample_name] = 0
		convert_sample_name = sample_name +"/"+ str(splicing_variant_name_dict[sample_name])
#		print("convert_sample_name: " , convert_sample_name, sep="\t")
		if line.startswith("CodingType"):
			continue
		else:
			print(convert_sample_name, line, sep="\t")

	
def novel_exon_convert2(f, f_2, f_3, f_4):

	def del_emp(l):
		l2 = []
		for i in range(len(l)):
			if not l[i] == "":
				l2.append(l[i])
		return l2

	gene_dict = {}
	f1 = open(f)
	for line in f1:
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
	f3 = open(f_3)
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
	f2 = open(f_2)
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

	f4 = open(f_4)
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

	
def novel_exon_convert3(f, f_2, int3, int4, f_5):

	def del_emp(l):
		l2 = []
		for i in range(len(l)):
			if not l[i] == "":
				l2.append(l[i])
		return l2

	gene_pos_dict = {}
	f1 = open(f)
	for line in f1:
		line = line.replace("\n","")
		line_l = line.split("\t")
		start_pos_l = del_emp(line_l[9].split(","))
		end_pos_l = del_emp(line_l[10].split(","))
		exon_num_l = line_l[-1].split(",")
		for i in range(len(start_pos_l)):
			k = line_l[12] +"/"+ exon_num_l[i]
			if not k in gene_pos_dict:
				gene_pos_dict[k] = [start_pos_l[i], end_pos_l[i]]
		
	novel_transcript_eva_dict = {}
	novel_transcript_convert_dict = {}
	f2 = open(f_2)
	for line in f2:
		line = line.replace("\n","")
		line_l = line.split("\t")
	#	print(line)
		if line.startswith("TranscriptID"):
			continue
		gene_info = line_l[0].split("/")
		novel_pos_l = line_l[2].split("/")
	#	print("novel_pos_l: ", novel_pos_l, sep="\t")
		novel_pos_l2 = novel_pos_l[1].split(",")
		exon_num_l = novel_pos_l[0].split(",")
		if "alt" in line_l[2] or "chrUn" in line_l[2] or "random" in line_l[2]:
	#		print("chrUn")
			novel_transcript_eva_dict[gene_info[2] +"/"+ line_l[1]] = False
			continue
	#	print("novel_pos_l2: ", novel_pos_l2, sep="\t")	
		start_l = []
		end_l = []
		novel_len_l = []
		for i in range(len(novel_pos_l2)):
	#		print(novel_pos_l2[i])	
			if "chr" in novel_pos_l2[i]:
				pos_l = novel_pos_l2[i].split("_")
	#			print("pos_l: ", pos_l, sep="\t")
				start_l.append(pos_l[1])
				end_l.append(pos_l[2])
				novel_len_l.append(int(pos_l[2]) - int(pos_l[1]) + 1)
			else:
				k = gene_info[2] + "/" + exon_num_l[i]
	#			print("gene_pos_dict[k]: ", gene_pos_dict[k], sep="\t")
				start_l.append(gene_pos_dict[k][0])
				end_l.append(gene_pos_dict[k][1])
	#	print("start_l: ", start_l, sep="\t")
	#	print("end_l: ", end_l, sep="\t")
	#	print("novel_len_l: ", novel_len_l, sep="\t")
		
		eva_novel_len = True
		for i in range(len(novel_len_l)):
			if int(novel_len_l[i]) < int(int3):
				eva_novel_len = False
	#	print("eva_novel_len: ", eva_novel_len, sep="\t")
		
		eva_intron_len = True
		for i in range(len(start_l)):
			if i == 0:
				continue
			intron_len = int(end_l[i]) - int(start_l[i-1]) + 1
	#		print("intron_len: ", intron_len, sep="\t")
			if intron_len < int(int4):
				eva_intron_len = False	
	#	print("eva_intron_len: ", eva_intron_len, sep="\t")
			
		if eva_novel_len == False or eva_intron_len == False:
	#		print("FilterOut")
			novel_transcript_eva_dict[gene_info[2] +"/"+ line_l[1]] = False
	#		print()
			continue
		novel_transcript_eva_dict[gene_info[2] +"/"+ line_l[1]] = True
		novel_transcript_convert_dict[gene_info[2] +"/"+ line_l[1]] = line_l[2]
	#	print()
		
	#for k,v in novel_transcript_eva_dict.items():
	#	print(k,v,sep="\t")#AFM/2,4,5,novel/s,k,k,chr4_73484454_73485677/*,*/*      True
		
	f5 = open(f_5)
	for line in f5:
		line = line.replace("\n","")
		line_l = line.split("\t")
	#	print(line)
		if line_l[2] == "UNANNOTATED":
	#		print(line)
	#		print(line_l[3] +"/"+ line_l[5])
	#		print(novel_transcript_eva_dict[line_l[3] +"/"+ line_l[5]])
			if novel_transcript_eva_dict[line_l[3] +"/"+ line_l[5]] == True:
	#			print(novel_transcript_convert_dict[line_l[3] +"/"+ line_l[5]])
				print("\t".join(line_l[0:5]), novel_transcript_convert_dict[line_l[3] +"/"+ line_l[5]], "\t".join(line_l[6:]), sep="\t")
	#		print()
	#	print()
		else:
			print(line)



