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
		start_pos_l = line_l[7].split(",")
		end_pos_l = line_l[8].split(",")
#		print("start_pos_l: ", start_pos_l, sep="\t")	
#		print("end_pos_l: ", end_pos_l, sep="\t")
		full_pos_l = []
		for i in range(len(start_pos_l)):
#			print(start_pos_l[i], end_pos_l[i], sep="\t")
			full_pos_l.append(start_pos_l[i])
			full_pos_l.append(end_pos_l[i])
#		print("full_pos_l: ", len(full_pos_l), full_pos_l, sep="\t")
		full_pos_s = ",".join(full_pos_l)
#		print("full_pos_s: ", full_pos_s, sep="\t")
		variant_dict[line_l[5]] = full_pos_s	
#		print("variant_dict[line_l[5]]: ", variant_dict[line_l[5]], sep="\t")
		print("\t".join(line_l[0:6]), "\t".join(line_l[9:]), sep="\t")
#		print()	
		continue
		
	else:
		start_pos_l = line_l[7].split(",")	
		end_pos_l = line_l[8].split(",")
#		print("start_pos_l: ", start_pos_l, sep="\t")
#		print("end_pos_l: ", end_pos_l, sep="\t")
		full_pos_l = []
		for i in range(len(start_pos_l)):
#			print(start_pos_l[i], end_pos_l[i], sep="\t")
			full_pos_l.append(start_pos_l[i])
			full_pos_l.append(end_pos_l[i])
#		print("full_pos_l: ", len(full_pos_l), full_pos_l, sep="\t")	
		full_pos_s = ",".join(full_pos_l)
#		print("full_pos_s: ", full_pos_s, sep="\t")	
		
		del full_pos_l[0]
		del full_pos_l[-1]
		full_pos_s2 = ",".join(full_pos_l)
#		print("full_pos_s2: ", full_pos_s2, sep="\t") 
			
		eva_overlap = False	
		for i in range(len(gene_variant_l)):
#			print(">>> gene_variant_l[i]: ", gene_variant_l[i], sep="\t")
			pre_gene_s = variant_dict[gene_variant_l[i]]
#			print("pre_gene_s: ", pre_gene_s, sep="\t")
			
			if eva_overlap == True:
				continue
	
			if full_pos_s2 in pre_gene_s:
#				print("OVERLAP", "MERGE")
				print("\t".join(line_l[0:4]), "PARTIAL", gene_variant_l[i], "\t".join(line_l[9:]), sep="\t")
				eva_overlap = True
#				print("--------------------")	
#			else:
#				print("NOT OVERLAP")		
#				print("--------------------")	
	
		if eva_overlap == False:
			gene_variant_l.append(line_l[5])
			variant_dict[line_l[5]] = full_pos_s
#			print("gene_variant_l: ", len(gene_variant_l), gene_variant_l, sep="\t")
#			print("variant_dict[line_l[5]]: ", variant_dict[line_l[5]], sep="\t")
			print("\t".join(line_l[0:6]), "\t".join(line_l[9:]), sep="\t")	
#		print()	 
		
#	print()
		
