import sys
import re

f = open(sys.argv[1])
for line in f:
	line = line.replace("\n","")
	line_l = line.split("\t")
#	print(line_l[12])
	if "/" in line_l[12]:
		gene_l = line_l[12].split("/")
		line_l[12] = gene_l[0]+gene_l[1]
		print("\t".join(line_l))
	else:
		print("\t".join(line_l))
