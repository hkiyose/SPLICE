import sys
import re

def get_report(f, f_2, f_3):

	num_line = 0
	num_reads = 0
	num_bases = 0
	f1 = open(f)
	for line in f1:
		num_line += 1
		line = line.replace("\n","")
		if num_line%4 == 2:
			num_reads += 1
			num_bases += len(line)

	num_line = 0
	num_filt_reads = 0
	num_filt_bases = 0
	longest_seq_num = 0
	f2 = open(f_2)
	for line in f2:
		num_line += 1
		line = line.replace("\n","")
		if num_line%4 == 2:
			num_filt_reads += 1
			num_filt_bases += len(line)
			if len(line) > longest_seq_num:
				longest_seq_num = len(line)

	mapped_read_num = 0
	mapped_base_num = 0
	total_num = 0
	M_num = 0
	R_num = 0        
	I_num = 0        
	D_num = 0
	#cigar_s = ""
	f3 = open(f_3)
	for line in f3:
		line = line.replace("\n","")
		line_l = line.split("\t")
	#	print(line)
		mapped_read_num += 1
		mapped_base_num += int(line_l[1])	
	#	print(line_l[10])
		cs_l = line_l[10].split("/")
		for k in range(len(cs_l)):	
			cs_l2 = cs_l[k].split('~')
	#		cigar_s = ""

			for i in range(len(cs_l2)):
				cs_l_type = re.findall('[:*+-]', cs_l2[i])
				cs_l_len = re.split('[:*+-]', cs_l2[i])
				if i == 0:
					del cs_l_type[0:2]
					del cs_l_len[0:3]
				else:
					del cs_l_len[0]
			
	#			print("cs_l_type: ", len(cs_l_type), cs_l_type )
	#			print("cs_l_len: ", len(cs_l_len), cs_l_len)
		
				for j in range(len(cs_l_type)):
					if cs_l_type[j] == ":":
	#					cigar_s += "M" * int(cs_l_len[j])
						M_num += int(cs_l_len[j])
					elif cs_l_type[j] == "*":
	#					cigar_s += "R" * len(cs_l_len[j])
						R_num += 1
					elif cs_l_type[j] == "+":
	#					cigar_s += "I" * len(cs_l_len[j])
						I_num += len(cs_l_len[j])
					else:#cs_l_type[j] == "-"
	#					cigar_s += "D" * len(cs_l_len[j])
						D_num += len(cs_l_len[j])
	total_num = M_num + R_num + I_num + D_num
	#print(cigar_s)
	#print(len(cigar_s))	
	#print()

	print("#total_reads", "total_bases", "filtered_reads", "filtered_bases", "longest_read_length", "mapped_reads", "mapped_bases", "average_length", "mismatch_rate%", "insertion_rate%", "deletino_rate%", sep="\t")
	print(num_reads, num_bases, num_filt_reads, num_filt_bases, longest_seq_num, mapped_read_num, mapped_base_num, round(mapped_base_num/mapped_read_num, 3), round(R_num/total_num*100,3), round(I_num/total_num*100,3), round(D_num/total_num*100,3), sep="\t")



