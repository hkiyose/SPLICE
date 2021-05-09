import sys
import re


f = open(sys.argv[1])
f2 = open(sys.argv[1])
gene_l = []
for line in f:
	line = line.replace("\n","")
	line_l = line.split("\t")
#	print(line)
	if not line_l[12] in gene_l and not "-" in line_l[12]:
		gene_l.append(line_l[12])
#		print(gene_l)

for line in f2:
	line = line.replace("\n","")
	line_l = line.split("\t")
	eva_cg = False	
	if "-" in line_l[12]:
#		print(">>> - in gene_name", line_l[12], line, sep="\t")
		
		gene_l2 = line_l[12].split("-")
		
		if gene_l2[0] == "AS":
			print(line)
			continue
		else:	
			for i in range(len(gene_l)):
				if line_l[12] == gene_l2[0] + "-" + gene_l[i]:
#					print(">>> conjoined_gene")
#					print(line)
					eva_cg = True
		
	if eva_cg == False:
		print(line)
