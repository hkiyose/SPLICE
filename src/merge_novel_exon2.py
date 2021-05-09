import sys
import re

novel_exon_dict = {}
f = open(sys.argv[1])
for line in f:
	line = line.replace("\n","")
	line_l = line.split("\t")
#	print(line)
	if line.startswith("Variant"):
		continue
	info_l = line_l[0].split(";")
#	print("info_l: ", info_l, sep="\t")
	chr_pos_l = info_l[3].split(",")
#	print("chr_pos_l: ", chr_pos_l, sep="\t")

	for i in range(len(chr_pos_l)):
		if not chr_pos_l[i] == "known" and not chr_pos_l[i] == "short" and not chr_pos_l[i] == "long" and not chr_pos_l[i] == "slong":
			chr_pos_l2 = chr_pos_l[i].split("_")
			if len(chr_pos_l2) == 5:
#				print("5!: ", line, sep="\t")
				k = chr_pos_l2[0]+"_"+chr_pos_l2[1]+"_"+chr_pos_l2[2], chr_pos_l2[3], chr_pos_l2[4]
				if k in novel_exon_dict:
					novel_exon_dict[k] += 1
				else:	
					novel_exon_dict[k] = 1
			elif len(chr_pos_l2) == 4:
#				print("4!: ", line, sep="\t")
				k = chr_pos_l2[0]+"_"+chr_pos_l2[1], chr_pos_l2[2], chr_pos_l2[3]
				if k in novel_exon_dict:
					novel_exon_dict[k] += 1
				else:
					novel_exon_dict[k] = 1
			else:#elif len(chr_pos_l2) == 3:
#				print("3!!!: ", line, sep="\t")
				k = chr_pos_l2[0], chr_pos_l2[1], chr_pos_l2[2]
				if k in novel_exon_dict:
					novel_exon_dict[k] += 1
				else:
					novel_exon_dict[k] = 1			
#	print()

for k,v in sorted(novel_exon_dict.items(), key=lambda x: (x[0][0], int(x[0][-2]), int(x[0][-1]))):
	print(k[0], k[1], k[2], v, sep="\t")
