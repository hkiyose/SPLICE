import sys
import re

merged_gene_bp_pos_dict = {}#merged_gene_bp_pos_dict[gene_A/gene_B] = [gene_A_bp_start_pos1, gene_A_bp_start_pos2, gene_B_bp_start_pos1, gene_B_bp_start_pos2]
f = open(sys.argv[1])
for line in f:
	line = line.replace("\n","")
	line_l = line.split("\t")
#	print(line)
	gene_l = line_l[2].split("/")
	bp_pos_l = line_l[4].split("/")
#	print(gene_l, bp_pos_l)
	
	k = "/".join(gene_l)
#	print("k: ", k, sep="\t")
	
	if k in merged_gene_bp_pos_dict:
#		print("<><><><><> indict <><><><><>")
#		print(merged_gene_bp_pos_dict[k])
#		print(bp_pos_l)
		for i in range(len(merged_gene_bp_pos_dict[k])):
			if int(merged_gene_bp_pos_dict[k][i][0]) - int(sys.argv[2]) <= int(bp_pos_l[0]) <= int(merged_gene_bp_pos_dict[k][i][1]) + int(sys.argv[2]) and int(merged_gene_bp_pos_dict[k][i][2]) - int(sys.argv[2]) <= int(bp_pos_l[1]) <= int(merged_gene_bp_pos_dict[k][i][3]) + int(sys.argv[2]) and line_l[3] == merged_gene_bp_pos_dict[k][i][7]:
#				print("overlap")
				if int(bp_pos_l[0]) < int(merged_gene_bp_pos_dict[k][i][0]):
					merged_gene_bp_pos_dict[k][i][0] = bp_pos_l[0]
				elif int(bp_pos_l[0]) > int(merged_gene_bp_pos_dict[k][i][1]):
					merged_gene_bp_pos_dict[k][i][1] = bp_pos_l[0]
		
				if int(bp_pos_l[1]) < int(merged_gene_bp_pos_dict[k][i][2]):
					merged_gene_bp_pos_dict[k][i][2] = bp_pos_l[1]
				elif int(bp_pos_l[1]) > int(merged_gene_bp_pos_dict[k][i][3]):
					merged_gene_bp_pos_dict[k][i][3] = bp_pos_l[1]
			
				merged_gene_bp_pos_dict[k][i][4] += int(line_l[0])
				merged_gene_bp_pos_dict[k][i][5] += float(line_l[1])
				merged_gene_bp_pos_dict[k][i][8] += "," + line_l[5]
			
#				print("merged", merged_gene_bp_pos_dict[k], sep="\t")
				break
			elif i == len(merged_gene_bp_pos_dict[k]) -1:
#				print("non_overlap")
				v = [bp_pos_l[0], bp_pos_l[0], bp_pos_l[1], bp_pos_l[1], int(line_l[0]), float(line_l[1]), line_l[2], line_l[3], line_l[5]]
	#			merged_gene_bp_pos_dict[k] = [bp_pos_l[0], bp_pos_l[0], bp_pos_l[1], bp_pos_l[1], int(line_l[0]), float(line_l[1]), line_l[2], line_l[3], line_l[5]]
				merged_gene_bp_pos_dict[k].append(v)
#				print(merged_gene_bp_pos_dict[k])
			
	else:
		v = [bp_pos_l[0], bp_pos_l[0], bp_pos_l[1], bp_pos_l[1], int(line_l[0]), float(line_l[1]), line_l[2], line_l[3], line_l[5]]
#		print("v: ", v, sep="\t")
		merged_gene_bp_pos_dict[k] = []
		merged_gene_bp_pos_dict[k].append(v)
#		print(merged_gene_bp_pos_dict[k])
	
#	print()

#for k,v in sorted(merged_gene_bp_pos_dict.items(), key=lambda x:[x[1][4],x[1][5]], reverse=True):
for k,v in merged_gene_bp_pos_dict.items():
#	print(k,v)
	for i in range(len(v)):
		if v[i][4] >= int(sys.argv[3]) and v[i][5] >= float(sys.argv[4]):
	 		print(v[i][4], round(v[i][5],3), k, v[i][7], v[i][0]+"-"+v[i][1]+"/"+v[i][2]+"-"+v[i][3], v[i][8], sep="\t")
