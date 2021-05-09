import sys
import re

f = open(sys.argv[1])
for line in f:
	line = line.replace("\n","")
	line_l = line.split("\t")
	if line.startswith("Variant"):
		continue
#	print(line)
	info_l = line_l[0].split(";")
#	print("info_l: ", info_l, sep="\t")
	
	annot_gene_l = info_l[1].split("/")
	annot_gene = "" 
	for i in range(len(annot_gene_l)):
#		print(annot_gene_l[i])
		if annot_gene_l[i] == "novel":
			continue
		else:
			if annot_gene == "":
				annot_gene = annot_gene_l[i]
			else:
				annot_gene += "/"+annot_gene_l[i]
	
	if annot_gene == "":
		annot_gene = "NOVEL"
#	print("annot_gene: ", annot_gene, sep="\t")
	
	eva_exon_len_l = info_l[3].split(",")
	convert_eva_exon_len_l = []
	for i in range(len(eva_exon_len_l)):
		if eva_exon_len_l[i] == "known":
			convert_eva_exon_len_l.append("k")
		elif eva_exon_len_l[i] == "short":
			convert_eva_exon_len_l.append("s")
		elif eva_exon_len_l[i] == "long":
			convert_eva_exon_len_l.append("l")
		elif eva_exon_len_l[i] == "slong":
			convert_eva_exon_len_l.append("sl")
		else:
			convert_eva_exon_len_l.append(eva_exon_len_l[i])

	exon_num_l = re.split('[,/]', info_l[2])
	convert_exon_num_l = []
#	print("exon_num_l: ", exon_num_l, sep="\t")
	exon_num_l2 = []
	for i in range(len(exon_num_l)):
		if exon_num_l[i] == "*":
			continue
		else:
			exon_num_l2.append(exon_num_l[i])
#	print("exon_num_l2: ", exon_num_l2, sep="\t")
	
	count = 0
	for i in range(len(convert_eva_exon_len_l)):
		if "_" in convert_eva_exon_len_l[i]:
			convert_exon_num_l.append("novel")		
		else:
			convert_exon_num_l.append(exon_num_l2[count])
			count += 1
			
#	print("convert_exon_num_l: ", convert_exon_num_l, sep="\t")
	
	print(annot_gene, "NOVEL", ",".join(convert_exon_num_l) +"/"+ ",".join(convert_eva_exon_len_l) +"/*,*/*", "\t".join(line_l[1:]), sep="\t")	
	
#	print()
