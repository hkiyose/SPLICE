import sys
import re

UTR3_merge_dict = {}

f = open(sys.argv[1])
for line in f:
	line = line.replace("\n","")
	line_l = line.split("\t")
#	print(line)
	if line.startswith("CodingType"):
		print(line)
		continue
	transcript_info = ""
	if "3UTR" in line_l[5]:
		transcript_info = ";".join(line_l[0:5])+";"+"-"
	else:
		transcript_info = ";".join(line_l[0:6])
#	print("transcript_info: " ,transcript_info, sep="\t")
	
	if not transcript_info in UTR3_merge_dict:
		UTR3_merge_dict[transcript_info] = []
		for reads in line_l[10:]:
			UTR3_merge_dict[transcript_info].append(int(reads))
	else:
#		print("MERGE!")
		for i in range(len(line_l[10:])):
			UTR3_merge_dict[transcript_info][i] += int(line_l[i+10])
#	print(UTR3_merge_dict[transcript_info])
#	print()

for k,v in UTR3_merge_dict.items():
	k_l = k.split(";")
	v_s = []
	for i in range(len(v)):
		v_s.append(str(v[i]))
	print("\t".join(k_l), "*", "*", "*", "*", "\t".join(v_s), sep="\t")
	
