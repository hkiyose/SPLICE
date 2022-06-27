import sys
import re

f = open(sys.argv[1])
gene_pos_dict = {}
for line in f:
	line = line.replace("\n", "")
	line_l = line.split("\t")
	if line.startswith("#"):
		continue
	pos_s = line_l[9] +"/"+ line_l[10]
#	print("pos_s: ", pos_s)
	if not line_l[12] in gene_pos_dict:
		gene_pos_dict[line_l[12]] = []
		gene_pos_dict[line_l[12]].append(pos_s)
		print(line)
	else:
		if not pos_s in gene_pos_dict[line_l[12]]:
			gene_pos_dict[line_l[12]].append(pos_s)
			print(line)
#		else:	
#			print(">>>", line)
