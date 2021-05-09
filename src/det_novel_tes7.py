import sys
import re

def del_terminal_index(l):
	if l[-1] == "":
		del l[-1]
	return l

f = open(sys.argv[1])
gene_exonNum_dict = {}
gene_strand_dict = {}
for line in f:
	line = line.replace("\n", "")
	line_l = line.split("\t")
	if not line_l[12] in gene_exonNum_dict:
		gene_exonNum_dict[line_l[12]] = []
		gene_exonNum_dict[line_l[12]].append(line_l[16])
		gene_strand_dict[line_l[12]] = line_l[3]
	else:
		gene_exonNum_dict[line_l[12]].append(line_l[16])

#for k,v in gene_exonNum_dict.items():
#	print(k,v,sep="\t")
#A1BG    ['3,7,9,11,14,16,17,19', '2,7,9,11,14,16,17,19', '10,
#A1BG-AS1        ['6,9,10,17', '4,8,10,16', '2,8,12', '3,9,13'

f3 = open(sys.argv[3])
annot_num_dict = {}
for line in f3:
	line = line.replace("\n", "")
	line_l = line.split("\t")
	annot_l = del_terminal_index(line_l[11].split("/"))#gene
	for i in range(len(annot_l)):
		if annot_l[i] in annot_num_dict:
			annot_num_dict[annot_l[i]] += 1
		else:
			annot_num_dict[annot_l[i]] = 1

f2 = open(sys.argv[2])
for line in f2:
	line = line.replace("\n", "")	
	line_l = line.split("\t")
#	print(line)
#	print("gene_strand_dict[line_l[2]]: ", gene_strand_dict[line_l[2]], sep="\t")
#	"+"：右のexonが3', "-"：左のexonが3'
	exonNum_l = line_l[6].split(",")
	eva_exon_len_l = del_terminal_index(line_l[7].split(","))
##	print("exonNum_l: ", exonNum_l, sep="\t")
##	print("eva_exon_len_l: ", eva_exon_len_l, sep="\t")

	if len(exonNum_l) == 1:#1exonのみ
#		print("continue")
#		print()
		continue
	
	if gene_strand_dict[line_l[2]] == "+":#Remove the exon on the 5' side. If there is still no correspondence, it can be judged as a new 3' side (polyA side).
		del exonNum_l[0]
		del eva_exon_len_l[0]
	else:
		del exonNum_l[-1]
		del eva_exon_len_l[-1]
	
##	print("masked_exonNum_l: ", exonNum_l, sep="\t")
##	print("masked_eva_exon_len_l: ", eva_exon_len_l, sep="\t")
	
	eva_match = False	
	for i in range(len(gene_exonNum_dict[line_l[2]])):
		ref_exonNum_l = gene_exonNum_dict[line_l[2]][i].split(",")
##		print("ref_exonNum_l: ", ref_exonNum_l, sep="\t")
		if set(exonNum_l) <= set(ref_exonNum_l):
			match_index_l = []
##			print(">>>ref_exonNum_l: ", ref_exonNum_l, sep="\t")
			eva_match = True
	
	if "short" in eva_exon_len_l or "long" in eva_exon_len_l or "slong" in eva_exon_len_l:
		eva_match = False

##	print("eva_match: ", eva_match, sep="\t")	
		
	if eva_match == False:
		print(line_l[1], round(int(line_l[1])/annot_num_dict[line_l[2]],3), line_l[2], line_l[3], line_l[4], line_l[5], line_l[6], line_l[7], line_l[8], line_l[9], sep="\t")	
##	print("")
#	else:#New candidate on the 5' side
		
		
