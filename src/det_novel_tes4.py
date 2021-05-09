import sys
import re

def del_terminal_index(l):
	if l[-1] == "":
		del l[-1]
	return l

f = open(sys.argv[1])

gene = ""
variant_geneNum_readLen_dict = {}#variant_geneNum_readLen_dict[A2M_9116157_9067707_*,7,10,11,*_*,known,known,known,*] = [1, a2550790-5011-4638-8078-951b86eb88f3]
for line in f:
	line = line.replace("\n", "")
	line_l = line.split("\t")
	
#	print(line)
	if gene == "":
		gene = line_l[2]
	if gene == line_l[2]:
		k = gene + ";" + line_l[3] +";"+ line_l[4] +";"+ line_l[5] +";"+ line_l[6] + ";" + line_l[7] +";"+ line_l[8]  
		if k in variant_geneNum_readLen_dict:
			variant_geneNum_readLen_dict[k][0] += 1
			variant_geneNum_readLen_dict[k][1] += ","+line_l[-1]
		else:
			variant_geneNum_readLen_dict[k] = [1, line_l[-1]]
	else:
		for k, v in variant_geneNum_readLen_dict.items():
			k_l = k.split(";")
			
			print(v[0], k_l[0], k_l[1], k_l[2], k_l[3], k_l[4], k_l[5], k_l[6], v[1], sep="\t")
		
		variant_geneNum_readLen_dict = {}
		
		gene = line_l[2]
		k = gene + ";" + line_l[3] +";"+ line_l[4] +";"+ line_l[5] +";"+ line_l[6] + ";" + line_l[7] +";"+ line_l[8]
		variant_geneNum_readLen_dict[k] = [1, line_l[-1]]

for k, v in variant_geneNum_readLen_dict.items():
	k_l = k.split(";")
	print(v[0], k_l[0], k_l[1], k_l[2], k_l[3], k_l[4], k_l[5], k_l[6], v[1], sep="\t")

