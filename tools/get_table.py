import sys
import re

def get_table(f, f_2, f_3, f_4, f_5):

	def get_eva_len_s(eva_len):
		eva_len_l = eva_len.split(",")
		eva_len_s = ""
		for i in range(len(eva_len_l)):
			if eva_len_l[i] == "*":
				eva_len_s += "*,"
			elif eva_len_l[i] == "known":
				eva_len_s += "k,"
			elif eva_len_l[i] == "short":
				eva_len_s += "s,"
			elif eva_len_l[i] == "long":
				eva_len_s += "l,"
			elif eva_len_l[i] == "slong":
				eva_len_s += "sl,"
		return eva_len_s[:-1]

	novel_read_dict = {}
	variant_expression_dict = {}
	variant_read_dict = {}
	variant_error_rate_dict = {}

	#variant_junction_match_rate_dict = {}
	f1 = open(f)
	for line in f1:
		line = line.replace("\n","")
		line_l = line.split("\t")
	#	print(line)
	#	if not int(line_l[0]) >= 2:
	#		continue
		eva_len_s = get_eva_len_s(line_l[8])	
	##	print("eva_len_s: ", eva_len_s, sep="\t")
		read_l = line_l[-1].split(",")

	#	if not int(line_l[10]) >= int(line_l[11]):
	#		continue

		for i in range(len(read_l)):
	#		print(read_l[i])
			novel_read_dict[read_l[i]] = 0
		k = line_l[2] +"/"+ "NOVEL" +"/"+ "combination" +"/"+ line_l[7] +"/"+ eva_len_s +"/"+ line_l[5] +","+ line_l[6] +"/"+ line_l[9] +"/"+ line_l[12]
	#CRP/NOVEL/combination/*,14,*/*,k,*/159712794,159714425/3,6,3,5/*,0,0,0,0,*	
		if not k in variant_expression_dict: 
			variant_expression_dict[k] = line_l[0]
			variant_read_dict[k] = line_l[-1]
	#		print(line_l[10], line_l[11])
			variant_error_rate_dict[k] = line_l[10]+","+line_l[11]
	##	else:	
	#		print(k, variant_expression_dict[k], variant_read_dict[k])
	##		print(">>>", line, sep="\t")
	#	print(k)
	#	print("variant_expression_dict: ", variant_expression_dict[k], sep="\t")
	#	print("variant_read_dict: ", variant_read_dict[k], sep="\t")
	#	print()

	f2 = open(f_2)
	for line in f2:
		line = line.replace("\n","")
		line_l = line.split("\t")
	##	print(line)
		eva_len_s = get_eva_len_s(line_l[8])
		read_l = line_l[-1].split(",")

	#	if not int(line_l[10]) >= int(line_l[11]):
	#		continue
		k = line_l[2] +"/"+ "NOVEL" +"/"+ "exon_length" +"/"+ line_l[7] +"/"+ eva_len_s +"/"+ line_l[5] +","+ line_l[6] +"/"+ line_l[9] +"/"+ line_l[12]
	#FTL/NOVEL/exon_length/*,3,4,*/*,k,s,*/48965609,48966583/3,7,1,1,4,6/*,0,0,0,3,0,0,*
		for i in range(len(read_l)):
			novel_read_dict[read_l[i]] = 0
		variant_expression_dict[k] = line_l[0]
		variant_read_dict[k] = line_l[-1]
		variant_error_rate_dict[k] = line_l[10]+","+line_l[11]
	##	print(line_l[2] +"/"+ "NOVEL" +"/"+ "exon_length" +"/"+ line_l[7] +"/"+ eva_len_s +"/"+ line_l[5] +","+ line_l[6] +"/"+ line_l[9]) 
	##	print("variant_expression_dict: ", variant_expression_dict[k], sep="\t")
	##	print("variant_read_dict: ", variant_read_dict[k], sep="\t")
	##	print()

	#B2M/NOVEL/exon_length/*,46,*/*,l,*/44711613,44717607/7,10,56,31

	f3 = open(f_3)
	for line in f3:
		line = line.replace("\n","")
		line_l = line.split("\t")
	##	print(line)
		eva_len_s = get_eva_len_s(line_l[7])
		read_l = line_l[-1].split(",")
		for i in range(len(read_l)):
			novel_read_dict[read_l[i]] = 0
		k = line_l[2] +"/"+ line_l[3] +"/"+ "3UTR" +"/"+ line_l[6] +"/"+ eva_len_s +"/"+ line_l[4] +","+ line_l[5] +"/*/" + line_l[8]
	#FTL/ENST00000331825.10,NM_000146.3/3UTR/2,3,4,5/k,k,k,s/48965609,48966583/*/0,0,0,0,0,0,0,-1
		variant_expression_dict[k] = line_l[0]
		variant_read_dict[k] = line_l[-1]
		variant_error_rate_dict[k] = "*,*"
	##	print("variant_expression_dict: ", variant_expression_dict[line_l[2] +"/"+ line_l[3] +"/"+ "3UTR" +"/"+ line_l[6] +"/"+ eva_len_s +"/"+ "*" +","+ "*" +"/"+"*"], sep="\t")
	##	print("variant_read_dict: ", variant_read_dict[line_l[2] +"/"+ line_l[3] +"/"+ "3UTR" +"/"+ line_l[6] +"/"+ eva_len_s +"/"+ "*" +","+ "*" +"/"+"*"], sep="\t")
	##	print()

	f4 = open(f_4)
	for line in f4:
		line = line.replace("\n","")
		line_l = line.split("\t")
	##	print(line_l)
		eva_len_s = get_eva_len_s(line_l[7])
		exon_num_l = line_l[6].split(",")
		
		read_l = line_l[-1].split(",")
		read_l2 = []
		for i in range(len(read_l)):
			if not read_l[i] in novel_read_dict:
				read_l2.append(read_l[i])

	##	print("len(read_l): ", len(read_l), sep="\t")
	##	print("len(read_l2): ", len(read_l2), sep="\t")

	#A1BG/ENST00000263100.7,ENST00000595014.1,ENST00000596924.1,ENST00000598345.1,NM_130786.3/KNOWN/5,6/s,s/58347029,58347353/*/*
	#	print(line_l[2] +"/"+ line_l[0] +"/"+ "" +"/"+ line_l[6] +"/"+ eva_len_s +"/"+ line_l[4] +","+ line_l[5] +"/"+ "*") 	
		
		if len(read_l2) == 0:
	#		print(line)
			continue
		
		for i in range(len(read_l2)):
			novel_read_dict[read_l2[i]] = 0

		if "KNOWN" in line_l[0]:
	##		print(line_l[2] +"/"+ line_l[3] +"/"+ "KNOWN" +"/"+ line_l[6] +"/"+ eva_len_s +"/"+ line_l[4] +","+ line_l[5] +"/"+ "*")
			variant_expression_dict[line_l[2] +"/"+ line_l[3] +"/"+ "KNOWN" +"/"+ line_l[6] +"/"+ eva_len_s +"/"+ line_l[4] +","+ line_l[5] +"/*/"+ line_l[8]] = len(read_l2)
			variant_read_dict[line_l[2] +"/"+ line_l[3] +"/"+ "KNOWN" +"/"+ line_l[6] +"/"+ eva_len_s +"/"+ line_l[4] +","+ line_l[5] +"/*/" + line_l[8]] = ",".join(read_l2)
			variant_error_rate_dict[line_l[2] +"/"+ line_l[3] +"/"+ "KNOWN" +"/"+ line_l[6] +"/"+ eva_len_s +"/"+ line_l[4] +","+ line_l[5] +"/*/" + line_l[8]] = "*,*"
	##		print("variant_expression_dict: ", variant_expression_dict[line_l[2] +"/"+ line_l[3] +"/"+ "KNOWN" +"/"+ line_l[6] +"/"+ eva_len_s +"/"+ line_l[4] +","+ line_l[5] +"/"+ "*"], sep="\t")
	##		print("variant_read_dict: ", variant_read_dict[line_l[2] +"/"+ line_l[3] +"/"+ "KNOWN" +"/"+ line_l[6] +"/"+ eva_len_s +"/"+ line_l[4] +","+ line_l[5] +"/"+ "*"], sep="\t")

		elif len(exon_num_l) == 1:
	##		print(line_l[2] +"/"+ line_l[3] +"/"+ "single_exon" +"/"+ line_l[6] +"/"+ eva_len_s +"/"+ line_l[4] +","+ line_l[5] +"/"+ "*")
			variant_expression_dict[line_l[2] +"/"+ line_l[3] +"/"+ "single_exon" +"/"+ line_l[6] +"/"+ eva_len_s +"/"+ line_l[4] +","+ line_l[5] +"/*/"+ line_l[8]] = len(read_l2)
			variant_read_dict[line_l[2] +"/"+ line_l[3] +"/"+ "single_exon" +"/"+ line_l[6] +"/"+ eva_len_s +"/"+ line_l[4] +","+ line_l[5] +"/*/"+ line_l[8]] = ",".join(read_l2)
			variant_error_rate_dict[line_l[2] +"/"+ line_l[3] +"/"+ "single_exon" +"/"+ line_l[6] +"/"+ eva_len_s +"/"+ line_l[4] +","+ line_l[5] +"/*/"+line_l[8]] = "*,*"
	##		print("variant_expression_dict: ", variant_expression_dict[line_l[2] +"/"+ line_l[3] +"/"+ "single_exon" +"/"+ line_l[6] +"/"+ eva_len_s +"/"+ line_l[4] +","+ line_l[5] +"/"+ "*"], sep="\t")
	##		print("variant_read_dict: ", variant_read_dict[line_l[2] +"/"+ line_l[3] +"/"+ "single_exon" +"/"+ line_l[6] +"/"+ eva_len_s +"/"+ line_l[4] +","+ line_l[5] +"/"+ "*"], sep="\t")

		else:
	##		print(line_l[2] +"/"+ line_l[3] +"/"+ "partially_known" +"/"+ line_l[6] +"/"+ eva_len_s +"/"+ line_l[4] +","+ line_l[5] +"/"+ "*")
			variant_expression_dict[line_l[2] +"/"+ line_l[3] +"/"+ "partially_known" +"/"+ line_l[6] +"/"+ eva_len_s +"/"+ line_l[4] +","+ line_l[5] +"/*/"+ line_l[8]] = len(read_l2)
			variant_read_dict[line_l[2] +"/"+ line_l[3] +"/"+ "partially_known" +"/"+ line_l[6] +"/"+ eva_len_s +"/"+ line_l[4] +","+ line_l[5] +"/*/"+line_l[8]] = ",".join(read_l2)
			variant_error_rate_dict[line_l[2] +"/"+ line_l[3] +"/"+ "partially_known" +"/"+ line_l[6] +"/"+ eva_len_s +"/"+ line_l[4] +","+ line_l[5] +"/*/"+ line_l[8]] = "*,*"
	##		print("variant_expression_dict: ", variant_expression_dict[line_l[2] +"/"+ line_l[3] +"/"+ "partially_known" +"/"+ line_l[6] +"/"+ eva_len_s +"/"+ line_l[4] +","+ line_l[5] +"/"+ "*"], sep="\t")
	##		print("variant_read_dict: ", variant_read_dict[line_l[2] +"/"+ line_l[3] +"/"+ "partially_known" +"/"+ line_l[6] +"/"+ eva_len_s +"/"+ line_l[4] +","+ line_l[5] +"/"+ "*"], sep="\t")
		
	##	print()

	gene_num_dict = {}
	f5 = open(f_5)
	for line in f5:
		line = line.replace("\n", "")
		line_l = line.split("\t")
		if line_l[11] == "Novel":
			continue	
		if line_l[11] in gene_num_dict:
			gene_num_dict[line_l[11]] += 1
		else:
			gene_num_dict[line_l[11]] = 1

	#read_count = 0
	for k,v in sorted(variant_expression_dict.items(), key=lambda x:-int(x[1])):
	#for k,v in variant_expression_dict.items():
		k_l = k.split("/")
	#	print(k,v,sep="\t")
	#	read_count += int(v)
		print(v, round(int(v)/gene_num_dict[k_l[0]]*100,3), k_l[0], k_l[1], k_l[2], k_l[3], k_l[4], k_l[5], k_l[6], variant_error_rate_dict[k], k_l[7], variant_read_dict[k], sep="\t")
	#	print(k_l[0], k_l[1], k_l[2], k_l[3], k_l[4], k_l[5], v, int(v)/gene_num_dict[k_l[0]]*100, variant_read_dict[k], sep="\t")
	#print(len(novel_read_dict))
	#print(read_count)
