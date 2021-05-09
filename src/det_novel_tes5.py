import sys
import re

def del_terminal_index(l):
	if l[-1] == "":
		del l[-1]
	return l

f = open(sys.argv[1])
gene_exonNum_dict = {}
gene_exonNum_variant_dict = {}
for line in f:
	line = line.replace("\n", "")
	line_l = line.split("\t")
	if not line_l[12] in gene_exonNum_dict:
		gene_exonNum_dict[line_l[12]] = []
		gene_exonNum_dict[line_l[12]].append(line_l[16])
	else:
		gene_exonNum_dict[line_l[12]].append(line_l[16])
	
	gene_exonNum_variant_dict[line_l[12] +"/"+ line_l[16]] = line_l[1] 

#for k,v in gene_exonNum_dict.items():
#	print(k,v,sep="\t")
#A1BG    ['3,7,9,11,14,16,17,19', '2,7,9,11,14,16,17,19', '10,
#A1BG-AS1        ['6,9,10,17', '4,8,10,16', '2,8,12', '3,9,13'

#1       A1BG    ENST00000598345.1       58346992        58346873        5       short,  69669822-56e5-4178-83ff-1f96e27b8c74
f2 = open(sys.argv[2])
for line in f2:
	line = line.replace("\n", "")	
	line_l = line.split("\t")
	eva_exon_len_l = del_terminal_index(line_l[6].split(","))
#	print(line)
#	print(eva_exon_len_l)
#	if line_l[2] == "NC":#181121確認
#		continue

	exonNum_l = line_l[5].split(",")
#	print("exonNum_l: ", exonNum_l)
#	print("eva_exon_len_l: ", eva_exon_len_l)

	match_variant = ""

	for i in range(len(gene_exonNum_dict[line_l[1]])):
		ref_exonNum_l = gene_exonNum_dict[line_l[1]][i].split(",")
#		print("ref_exonNum_l: ", ",".join(ref_exonNum_l), "variant: ", gene_exonNum_variant_dict[line_l[1] +"/"+ ",".join(ref_exonNum_l)] , sep="\t")
		if line_l[5] == ",".join(ref_exonNum_l):
			if not "short" in line_l[6] and not "long" in line_l[6]:
#				print(">>>PERFECT MATCH")
				match_variant = gene_exonNum_variant_dict[line_l[1] +"/"+ ",".join(ref_exonNum_l)]
				break
				

#	print("match_variant: ", match_variant, sep="\t")
	
	if match_variant == "":
		print("NOVEL: ", line, sep="\t")
	else:	
		print("KNOWN: ", line_l[0], line_l[1], match_variant, line_l[3], line_l[4], line_l[5], line_l[6], line_l[7], line_l[8], sep="\t")

#	print("")	

