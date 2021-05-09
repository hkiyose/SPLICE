import sys
import re

def del_terminal_index(l):
	if l[-1] == "":
		del l[-1]
	return l

f = open(sys.argv[1])

gene = ""
gene_num_dict = {}
variant_geneNum_readLen_dict = {}#variant_geneNum_readLen_dict[A2M_9116157_9067707_*,7,10,11,*_*,known,known,known,*] = [1, a2550790-5011-4638-8078-951b86eb88f3]
junction_match_num_dict = {}
genome_transcriptome_match_num_dict = {}
for line in f:
	line = line.replace("\n", "")
	line_l = line.split("\t")
	
#	print(line)
	
	if line_l[2] in gene_num_dict:
		gene_num_dict[line_l[2]] += 1
	else:
		gene_num_dict[line_l[2]] = 1
		
	junction_match_count_l = line_l[9].split(",")

	if gene == "":
		gene = line_l[2]
	
	gap_l2 = []
	gap_l = line_l[-1].split(",")
	for i in range(len(gap_l)):
		if i == 0 or i == len(gap_l)-1:
			gap_l2.append("*")	
			continue
		else:
			gap_l2.append(gap_l[i])
		
#	print(",".join(gap_l2))
	
	if not line_l[11] == "*":
		if int(line_l[10]) < int(line_l[11]):
			continue

	if gene == line_l[2]:
	
		k = gene + ";" + "*" +";"+ line_l[4] +";"+ line_l[5] +";"+ line_l[6] + ";" + line_l[7] + ";" + line_l[8] + ";" + ",".join(gap_l2)
		
		if k in variant_geneNum_readLen_dict:
			variant_geneNum_readLen_dict[k][0] += 1
			variant_geneNum_readLen_dict[k][1] += ","+line_l[12]
	
			if not junction_match_count_l[0] == "*":
				for i in range(len(junction_match_count_l)):
					junction_match_num_dict[k][i] += int(junction_match_count_l[i])
			
			genome_transcriptome_match_num_dict[k][0] += int(line_l[10])
			if not line_l[11] == "*":
				genome_transcriptome_match_num_dict[k][1] += int(line_l[11])
			else:
				genome_transcriptome_match_num_dict[k][1] += 0
			
		else:
			variant_geneNum_readLen_dict[k] = [1, line_l[12]]
			
			if not junction_match_count_l[0] == "*":
				junction_match_num_dict[k] = [0]*len(junction_match_count_l)
				for i in range(len(junction_match_count_l)):
					junction_match_num_dict[k][i] += int(junction_match_count_l[i])
			else:
				junction_match_num_dict[k] = "*"
		
			if not line_l[11] == "*":
				genome_transcriptome_match_num_dict[k] = [int(line_l[10]), int(line_l[11])]
			else:
				genome_transcriptome_match_num_dict[k] = [int(line_l[10]), 0]
				
	
	else:
		for k, v in variant_geneNum_readLen_dict.items():
			k_l = k.split(";")
			junction_match_num_l2 = []
			if junction_match_num_dict[k][0] == "*":
				junction_match_num_l2.append("*")
			else:
				for i in range(len(junction_match_num_dict[k])):
					junction_match_num_l2.append(str(100-int(junction_match_num_dict[k][i]/int(v[0])/5*100)))
			
			genome_transcriptome_match_num_l = ["",""]
			genome_transcriptome_match_num_l[0] = str(int(int(genome_transcriptome_match_num_dict[k][0])/int(v[0])))
			genome_transcriptome_match_num_l[1] = str(int(int(genome_transcriptome_match_num_dict[k][1])/int(v[0])))
			
			print(v[0], round(v[0]/gene_num_dict[gene],2), k_l[0], k_l[1], k_l[2], k_l[3], k_l[4], k_l[5], k_l[6], ",".join(junction_match_num_l2), genome_transcriptome_match_num_l[0], genome_transcriptome_match_num_l[1], k_l[7], v[1], sep="\t")
		
		variant_geneNum_readLen_dict = {}
		
		gene = line_l[2]
		k = gene + ";" + "*" +";"+ line_l[4] +";"+ line_l[5] +";"+ line_l[6] + ";" + line_l[7] + ";" + line_l[8] +";"+ ",".join(gap_l2)
		variant_geneNum_readLen_dict[k] = [1, line_l[12]]
			
		if not junction_match_count_l[0] == "*":
			junction_match_num_dict[k] = [0]*len(junction_match_count_l)
			for i in range(len(junction_match_count_l)):
				junction_match_num_dict[k][i] += int(junction_match_count_l[i])
		else:
			junction_match_num_dict[k] = "*"
			
		if not line_l[11] == "*":
			genome_transcriptome_match_num_dict[k] = [int(line_l[10]), int(line_l[11])]
		else:
			genome_transcriptome_match_num_dict[k] = [int(line_l[10]), 0]

for k, v in variant_geneNum_readLen_dict.items():
	k_l = k.split(";")
	junction_match_num_l2 = []
	if junction_match_num_dict[k][0] == "*":
		junction_match_num_l2.append("*")
	else:
		for i in range(len(junction_match_num_dict[k])):
			junction_match_num_l2.append(str(100-int(junction_match_num_dict[k][i]/int(v[0])/5*100)))
	
	genome_transcriptome_match_num_l = ["",""]
	genome_transcriptome_match_num_l[0] = str(int(int(genome_transcriptome_match_num_dict[k][0])/int(v[0])))
	genome_transcriptome_match_num_l[1] = str(int(int(genome_transcriptome_match_num_dict[k][1])/int(v[0])))
		
	print(v[0], round(v[0]/gene_num_dict[gene],2), k_l[0], k_l[1], k_l[2], k_l[3], k_l[4], k_l[5], k_l[6], ",".join(junction_match_num_l2), genome_transcriptome_match_num_l[0], genome_transcriptome_match_num_l[1], k_l[7], v[1], sep="\t")

		
