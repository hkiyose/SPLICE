import sys
import re

f = open(sys.argv[1])
pre_read_name = ""
pre_match_count = 0
for line in f:
	line = line.replace("\n","")
	line_l = line.split("\t")
#	print(line)	
	if line.startswith("@"):
		continue
	if line_l[2] == "*":
		continue
	
#	print(line)
		
	for i in range(len(line_l)):
		if line_l[i].startswith("cs:"):
			cs = line_l[i]

	cs_l_type = re.findall('[:*+-]', cs)
	cs_l_len = re.split('[:*+-]', cs)
	del cs_l_type[0:2]
	del cs_l_len[0:3]
		
#	print("cs_l_type: ", len(cs_l_type), cs_l_type, sep="\t")
#	print("cs_l_len: ", len(cs_l_len), cs_l_len, sep="\t")

	cigar_s = ""
	for j in range(len(cs_l_type)):
		if cs_l_type[j] == ":":
			cigar_s += "M" * int(cs_l_len[j])
		elif cs_l_type[j] == "*":	
			cigar_s += "R"
		elif cs_l_type[j] == "+":
			cigar_s += "I" * len(cs_l_len[j])
		else:#cs_l_type[j] == "-"
			cigar_s += "D" * len(cs_l_len[j])
	
	if pre_read_name == "":
		pre_read_name = line_l[1]
		pre_match_count = cigar_s.count("M")		
		print(cigar_s.count("M"), line, sep="\t")
	elif pre_read_name == line_l[1]:
#		if pre_match_count == cigar_s.count("M"):
		print(cigar_s.count("M"), line, sep="\t")
#		else:
#			continue
	else:
		pre_read_name = line_l[1]
		pre_match_count = cigar_s.count("M")
		print(cigar_s.count("M"), line, sep="\t")
	
#	print()		
