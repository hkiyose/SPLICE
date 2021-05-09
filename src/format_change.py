import sys
import re

def del_terminal_index(l):
	if l[-1] == "":
		del l[-1]
	return l

f = open(sys.argv[1])
false_positive_num = 0
gene_num_dict = {}
novel_variant_dict = {}#dict[[gene/exon_num/eva_exon_len]] = read_num
read_id_dict = {}#dict[[gene/exon_num/eva_exon_len]] = read_ID
variant_name_dict = {}#dict[[gene/exon_num/eva_exon_len]] = variant_name
junction_match_count_dict = {}
genome_transcriptome_match_num_dict = {}
for line in f:
	line = line.replace("\n", "")
	line_l = line.split("\t")
##	print(line)
	annot_gene_l = del_terminal_index(line_l[11].split("/"))
	exon_num_l = del_terminal_index(line_l[12].split("/"))
	eva_exon_len_l = del_terminal_index(line_l[13].split("/"))
	tsstes_l = del_terminal_index(line_l[15].split("/"))
	variant_l = del_terminal_index(line_l[14].split("/"))
	junction_match_count_l = line_l[16].split(",") 
	
##	print("annot_gene_l: ", annot_gene_l)
#	print("exon_num_l: ", exon_num_l)
#	print("eva_exon_len_l: ", eva_exon_len_l)
##	print("tsstes_l: ", tsstes_l)		
#	print("variant_l: ", variant_l)	
##	print("junction_match_count_l: ", junction_match_count_l)
#	print("")
	
	annot_gene_num = 0
	for i in range(len(annot_gene_l)):
		if annot_gene_l[i] == "novel":
			annot_gene_num += 1
#			print(line)
			continue
		else:
			annot_gene_num += 1
			if not annot_gene_l[i] in gene_num_dict:
				gene_num_dict[annot_gene_l[i]] = 1
			else:
				gene_num_dict[annot_gene_l[i]] += 1
	
##	if annot_gene_num >= 2:
##		print("annot_gene_num >= 2: ", line)
##		print()
##		continue
	
	if annot_gene_num == 1:#fusionを除く
		for i in range(len(annot_gene_l)):
			if annot_gene_l[i] == "FP1" or annot_gene_l[i] == "FP2" or annot_gene_l[i] == "RP1" or annot_gene_l[i] == "RP2" or annot_gene_l[i] == "novel" or exon_num_l == []:
##				print("NOVEL: ", line)
				continue
			else:
				key_s = annot_gene_l[i] +"/"+ exon_num_l[i] +"/"+ eva_exon_len_l[i]#novel_variant_dictのkey
#				print(line)
#				print("key_s: ", key_s)
				
				tsstes_l_l = del_terminal_index(tsstes_l[i].split(","))
				variant_l_l = del_terminal_index(variant_l[i].split(","))
##				print("tsstes_l_l", tsstes_l_l)
##				print("varint_l_l", variant_l_l)
#				print("")
				
				variant2_l = [] 
				
				if "TSS+TES" in tsstes_l[i] or "TSS" in tsstes_l[i] or "TES" in tsstes_l[i]:
					key_s += "/"+"*"+"/"+",".join(variant_l_l)
				else:
					key_s += "/"+"NOVEL"+"/"+"*"
				
##				print("key_s: ", key_s)
				if not key_s in novel_variant_dict:
					novel_variant_dict[key_s] = 1
					read_id_dict[key_s] = line_l[0]
					if not junction_match_count_l[0] == "*":
						junction_match_count_dict[key_s] = [0]*len(junction_match_count_l)
					else:
						junction_match_count_dict[key_s] = "*"
					if line_l[18] == "*":
						genome_transcriptome_match_num_dict[key_s] = [int(line_l[17]), 0]
					else:
						genome_transcriptome_match_num_dict[key_s] = [int(line_l[17]), int(line_l[18])]
				else:
					novel_variant_dict[key_s] += 1
					read_id_dict[key_s] += ","+line_l[0]
					genome_transcriptome_match_num_dict[key_s][0] += int(line_l[17])
					if line_l[18] =="*":
						genome_transcriptome_match_num_dict[key_s][1] += 0
					else:
						genome_transcriptome_match_num_dict[key_s][1] += int(line_l[18])
				
				if not junction_match_count_l[0] == "*":
					for i in range(len(junction_match_count_l)):
						junction_match_count_dict[key_s][i] += int(junction_match_count_l[i])	
					
				
#				print("junction_match_count_dict[key_s]: ", junction_match_count_dict[key_s])		
#				print("genome_transcriptome_match_num_dict[key_s]: ", genome_transcriptome_match_num_dict[key_s])				

##				print("novel_variant_dict[key_s]: ", novel_variant_dict[key_s])
##				print("read_id_dict[key_s]: ", read_id_dict[key_s])
#		print(line)
#	else:
##		print()
	
##	print("")

#annot_gene_l:  ['HP']
#tsstes_l:  ['*']
#junction_match_count_l:  ['5', '5', '5', '5', '5', '5', '5', '5', '5', '3', '5', '5']
#0004419d-78bd-4f3b-9cd5-aa4def155248    1494    16      chr16   55,984,1060,1162,1237,1339,1421 983,1059,1161,1236,1338,1420,1452       72054626,72056161,72056530,72057392,72058254,72059114,72060112  72054657,72056243,72056631,72057466,7205835
#tsstes_l_l ['*']
#varint_l_l ['Novel']
#key_s:  HP/5,9,11,13,16,19,27/known,known,known,known,known,known,known,/NOVEL/*
#junction_match_count_dict[key_s]:  [25, 25, 25, 25, 25, 24, 24, 23, 25, 22, 24, 21]
#novel_variant_dict[key_s]:  5
#read_id_dict[key_s]:  00005722-b5c1-4fe1-9daa-fdd321f592e2,000398ba-2e86-42b8-8f3b-5301519ce8f8,0003a5c8-cec5-4b1f-80fe-b32ea66edc4b,00043315-80fa-4282-8884-1b383831a057,0004419d-78bd-4f3b-9cd5-aa4def155248

for k,v in sorted(novel_variant_dict.items(), key=lambda x: -x[1]):
	key_l = k.split("/")
#	print(key_l, v, sep="\t")
#	print("junction_match_count_dict[k]: ", junction_match_count_dict[k])
#	print("genome_transcriptome_match_num_dict[k]: ", genome_transcriptome_match_num_dict[k])

	junction_match_count_l2 = []
	if junction_match_count_dict[k][0] == "*":
		junction_match_count_l2.append("*") 
	else:
		for i in range(len(junction_match_count_dict[k])):
			junction_match_count_l2.append(str(100-int(junction_match_count_dict[k][i]/int(v)/5*100)))
		
	genome_transcriptome_match_num_l = ["",""]
	genome_transcriptome_match_num_l[0] = str(int(int(genome_transcriptome_match_num_dict[k][0])/int(v)))
	genome_transcriptome_match_num_l[1] = str(int(int(genome_transcriptome_match_num_dict[k][1])/int(v)))
		
#	print(junction_match_count_l2)
#	print(genome_transcriptome_match_num_l)
		
		
	eva_exon_len_l2 = del_terminal_index(key_l[2].split(","))
#	print("eva_exon_len_l2: ", eva_exon_len_l2)
	eva_exon_len_middle_slice_l = eva_exon_len_l2[1:-1]
	if not key_l[3] == "NOVEL":
		
		if not "short" in key_l[2] and not "long" in key_l[2] and not "slong" in key_l[2]:
			print(v, round(float(int(v)/int(gene_num_dict[key_l[0]])),2), key_l[0],key_l[4], "KNOWN", key_l[1], key_l[2], ",".join(junction_match_count_l2), ",".join(genome_transcriptome_match_num_l), read_id_dict[k], sep="\t")

		else:
			print(v, round(float(int(v)/int(gene_num_dict[key_l[0]])),2), key_l[0],key_l[4], "NOVEL", key_l[1], key_l[2], ",".join(junction_match_count_l2), ",".join(genome_transcriptome_match_num_l), read_id_dict[k], sep="\t")
	else:
		print(v, round(float(int(v)/int(gene_num_dict[key_l[0]])),2), key_l[0],key_l[4], "NOVEL", key_l[1], key_l[2], ",".join(junction_match_count_l2), ",".join(genome_transcriptome_match_num_l), read_id_dict[k], sep="\t")

