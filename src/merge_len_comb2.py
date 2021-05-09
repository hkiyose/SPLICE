import sys
import re

pre_gene = ""
f = open(sys.argv[1])
for line in f:
	line = line.replace("\n","")
	line_l = line.split("\t")
#	print(line)
	if line_l[3] != pre_gene:
		pre_gene = line_l[3]
		gene_variant_l =[]
		variant_dict = {}
		gene_variant_l.append(line_l[5])
		pos_s = line_l[7]
#		print("pos_s: ", pos_s, sep="\t")
		variant_dict[line_l[5]] = pos_s
#		print("variant_dict[line_l[5]]: ", variant_dict[line_l[5]], sep="\t")
		print("\t".join(line_l[0:6]), "\t".join(line_l[8:]), sep="\t")
#		print()
		continue
	else:
		pos_s = line_l[7]
#		print("pos_s: ", pos_s, sep="\t")
		eva_overlap = False
		for i in range(len(gene_variant_l)):
#			print(">>> gene_variant_l[i]: ", gene_variant_l[i], sep="\t")
			pre_gene_s = variant_dict[gene_variant_l[i]]
#			print("pre_gene_s: ", pre_gene_s, sep="\t")
			
			if eva_overlap == True:
				continue
			
			if pos_s in pre_gene_s:
#				print("OVERLAP", "MERGE")
				print("\t".join(line_l[0:4]), "PARTIAL", gene_variant_l[i], "\t".join(line_l[8:]), sep="\t")	
				eva_overlap = True
#				print("--------------------")
#			else:
#				print("NOT OVERLAP")
#				print("--------------------")
		if eva_overlap == False:
			gene_variant_l.append(line_l[5])
			variant_dict[line_l[5]] = pos_s
#			print("gene_variant_l: ", len(gene_variant_l), gene_variant_l, sep="\t")
#			print("variant_dict[line_l[5]]: ", variant_dict[line_l[5]], sep="\t")
			print("\t".join(line_l[0:6]), "\t".join(line_l[8:]), sep="\t")
#		print()
#	print()
