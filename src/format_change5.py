import sys
import re


protein_coding_dict = {}
f = open(sys.argv[1])
for line in f:
	line = line.replace("\n","")
	line_l = line.split("\t")
#       print(line)
	if line_l[1] == "protein_coding":
		protein_coding_dict[line_l[0]] = 0

f2 = open(sys.argv[2])
for line in f2:
	line = line.replace("\n","")
	line_l = line.split("\t")
	if not line_l[0] in protein_coding_dict and line_l[1] == "protein_coding":
		protein_coding_dict[line_l[0]] = 0

gene_chr_dict = {}
f3 = open(sys.argv[3])
for line in f3:
	line = line.replace("\n", "")
	line_l = line.split("\t")
#	print(line)
	if not line_l[12] in gene_chr_dict:
		gene_chr_dict[line_l[12]] = line_l[2]
	
	if line_l[2] == "chrM":
		if line_l[13] == "cmpl" and line_l[14] == "cmpl":
			protein_coding_dict[line_l[12]] = 0

f4 = open(sys.argv[4])
for line in f4:
	line = line.replace("\n","")
	line_l = line.split("\t")	
#	print(line)
	
	if line_l[0] == "CodingType":
#		print(line)
		print("CodingType", "EvaKnown", "EvaFullLength", "Chr", "Gene", "SplicingVariant", "NovelSplicingVariantInfo", "\t".join(line_l[10:]), sep="\t")
		continue
	
	if line_l[3] in protein_coding_dict:
		coding_type = "CODING"
	else:
		coding_type = "NON-CODING"
	
#	print("coding_type: ", coding_type, sep="\t")
	
#	print("gene_chr_dict[line_l[3]]: ", gene_chr_dict[line_l[3]], "gene_strand_dict[line_l[3]]: ", gene_strand_dict[line_l[3]], sep="\t")	
	
	if line_l[1] != "mtDNA":
		print(coding_type, line_l[1], line_l[2],  gene_chr_dict[line_l[3]], "\t".join(line_l[3:6]), "\t".join(line_l[10:]), sep="\t")	
	else:
		print(coding_type, "KNOWN", line_l[2], gene_chr_dict[line_l[3]], "\t".join(line_l[3:6]), "\t".join(line_l[10:]), sep="\t")
	
#	print()

