import sys
import re

f = open(sys.argv[1])
for line in f:
	line = line.replace("\n", "")
	line_l = line.split("\t")
#	print(line)
	if "/" in line_l[2]:
		continue	
	for i in range(len(line_l)):
		if line_l[i].startswith("cs:"):
			cs = line_l[i]
#	print(cs)

	cs_l = cs.split('~')
	cigar_s = ""
		
	for i in range(len(cs_l)):
		cs_l_type = re.findall('[:*+-]', cs_l[i])
		cs_l_len = re.split('[:*+-]', cs_l[i])	
		if i == 0:
			del cs_l_type[0:2]
			del cs_l_len[0:3]
		else:
			del cs_l_len[0]
		
#		print("cs_l_type: ", len(cs_l_type), cs_l_type )
#		print("cs_l_len: ", len(cs_l_len), cs_l_len)

		for j in range(len(cs_l_type)):
			if cs_l_type[j] == ":":
				cigar_s += "M" * int(cs_l_len[j])
			elif cs_l_type[j] == "*":
				cigar_s += "R"
			elif cs_l_type[j] == "+":
				cigar_s += "I" * len(cs_l_len[j])
			else:#cs_l_type[j] == "-"
				cigar_s += "D" * len(cs_l_len[j])
#	print(cigar_s)	
	print(cigar_s.count("M"),cigar_s.count("M")/int(line_l[1]), line, cigar_s, sep="\t")	
#	print("")
