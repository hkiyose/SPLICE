import sys
import re

f = open(sys.argv[1])
for line in f:
	line = line.replace("\n","")
	line_l = line.split("\t")
#	print(line)
	if "/" in line_l[2]:
		continue
##	print(line)
	cs_l = line_l[10].split('~')
	cigar_s = ""
	read_strat_pos_l = line_l[4].split(",")
	read_end_pos_l = line_l[5].split(",")
	read_mapped_len = (int(line_l[1]) - (int(read_strat_pos_l[0]) -1) - (int(line_l[1]) - int(read_end_pos_l[-1])))
##	print(read_mapped_len)
	
	for i in range(len(cs_l)):
		cs_l_type = re.findall('[:*+-]', cs_l[i])
		cs_l_len = re.split('[:*+-]', cs_l[i])
		if i == 0:
			del cs_l_type[0:2]
			del cs_l_len[0:3]
#		else:
#			del cs_l_len[0]
		
#		print("cs_l_type: ", len(cs_l_type), cs_l_type, sep="\t")
#		print("cs_l_len: ", len(cs_l_len), cs_l_len, sep="\t")

#		cigar_s = ""
		for j in range(len(cs_l_type)):
			if len(cs_l_type) == len(cs_l_len):
				if cs_l_type[j] == ":":
					cigar_s += "M" * int(cs_l_len[j])
				elif cs_l_type[j] == "*":	
					cigar_s += "R"
				elif cs_l_type[j] == "+":
					cigar_s += "I" * len(cs_l_len[j])
				else:#cs_l_type[j] == "-"
					cigar_s += "D" * len(cs_l_len[j])
			else:
				if j == 0:
#					cigar_s += ";" + cs_l_len[j] + ";"
					cigar_s += ";"
					if cs_l_type[j] == ":":
						cigar_s += "M" * int(cs_l_len[j+1])
					elif cs_l_type[j] == "*":
						cigar_s += "R"
					elif cs_l_type[j] == "+":
						cigar_s += "I" * len(cs_l_len[j+1])
					else:
						cigar_s += "D" * len(cs_l_len[j+1])
					
				else:	
					if cs_l_type[j] == ":":
						cigar_s += "M" * int(cs_l_len[j+1]) 
					elif cs_l_type[j] == "*":
						cigar_s += "R"
					elif cs_l_type[j] == "+":
						cigar_s += "I" * len(cs_l_len[j+1])
					else:
						cigar_s += "D" * len(cs_l_len[j+1])
#		print(cigar_s)	
#		print(cigar_s.count("M") + cigar_s.count("I") + cigar_s.count("R"))
##	print(cigar_s)
##	print("M+I+R: ", cigar_s.count("M") + cigar_s.count("I") + cigar_s.count("R"), sep="\t")
##	print(float(cigar_s.count("M")/(read_mapped_len + cigar_s.count("D"))*100))
	cigar_l = cigar_s.split(";")
	junction_match_num_l = []
	for i in range(len(cigar_l)):
#		print(">>", i, cigar_l[i])
		if i == 0:
			continue
##			print(cigar_l[i][0:5])
		else:
##			print(cigar_l[i-1][-5::], cigar_l[i][0:5])
			junction_match_num_l.append(str(cigar_l[i-1][-5::].count("M")))
			junction_match_num_l.append(str(cigar_l[i][0:5].count("M")))
#			junction_seq = cigar_l[i-1][-5::]+cigar_l[i][0:5]
#			print(junction_seq)
	if len(junction_match_num_l) == 0:
		print(line, "*", cigar_s.count("M"), sep="\t")
#		print("*")
	else:
		print(line, ",".join(junction_match_num_l), cigar_s.count("M"), sep="\t")
#		print(",".join(junction_match_num_l))
##	print()
