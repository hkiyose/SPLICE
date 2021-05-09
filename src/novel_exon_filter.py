import sys
import re

protein_coding_dict = {}
f = open(sys.argv[1])
for line in f:
	line = line.replace("\n","")
	line_l = line.split("\t")
#	print(line)
	if line_l[1] == "protein_coding":
		protein_coding_dict[line_l[0]] = 0

f2 = open(sys.argv[2])
for line in f2:
	line = line.replace("\n","")
	line_l = line.split("\t")
	if not line_l[0] in protein_coding_dict and line_l[1] == "protein_coding":
		protein_coding_dict[line_l[0]] = 0

#for k,v in protein_coding_dict.items():
#	print(k,v,sep="\t")

f3 = open(sys.argv[3])
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









