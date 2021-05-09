import sys
import re

def del_terminal_index(l):
	if l[-1] == "":
		del l[-1]
	return l

f = open(sys.argv[1])
convert_gene_into_variant_dict = {}
variant_strand_direction_dict = {}
isoNum_dict = {}
isoNum_uniq_dict = {}

for line in f:
	line = line.replace("\n","")
	line_l = line.split("\t")
	if not len(line_l) >= 2:
		continue 
	convert_gene_into_variant_dict[line_l[12],line_l[16]] = line_l[1]
	variant_strand_direction_dict[line_l[1]] = line_l[3]
    
	if not line_l[12] in isoNum_dict:
		exon_num_list = [line_l[16]]
		isoNum_dict[str(line_l[12])] = exon_num_list
	else:
		exon_num_list = [line_l[16]]
		isoNum_dict[str(line_l[12])] += exon_num_list

for k ,v in isoNum_dict.items():
	uniq_v = list(set(v))
	isoNum_uniq_dict[k] = uniq_v

isoNum_uniq_dict['novel'] = '*'
isoNum_uniq_dict['FP1'] = '1'
isoNum_uniq_dict['RP1'] = '1'
isoNum_uniq_dict['FP2'] = '1'
isoNum_uniq_dict['RP2'] = '1'
isoNum_uniq_dict['P1'] = '1'
isoNum_uniq_dict['P2'] = '1'
#print(isoNum_uniq_dict['VEGFA'])
#print(isoNum_dict['VEGFA'])

f2 = open(sys.argv[2])

for line in f2:
	line = line.replace("\n","")
	line_l = line.split("\t")
   
#	print(line) 
 
	gene_l = line_l[11].split("/")
#	print(gene_l)
	if len(gene_l) > 1:
		del gene_l[-1]
   
	exonNum_l = line_l[12].split("/")
#	print(exonNum_l)
	if len(exonNum_l) > 1:
		del exonNum_l[-1]

	ref_exonNum_l = []
	for i in range(len(gene_l)):
		if gene_l[i] == "TRAV14":#refの遺伝子名に"/"が含まれる
			ref_exonNum_l.append(isoNum_uniq_dict[gene_l[i]+"/DV4"]) 
		elif gene_l[i] == "THRA1":
			ref_exonNum_l.append(isoNum_uniq_dict[gene_l[i]+"/BTR"])
		elif gene_l[i] in isoNum_uniq_dict:
			ref_exonNum_l.append(isoNum_uniq_dict[gene_l[i]])
	
#	print("gene_l: ", gene_l)
#	print("exonNum_l: ", exonNum_l)
#	print("ref_exonNum_l: ", ref_exonNum_l)

#	if line_l[9] in isoNum_uniq_dict[line_l[8]]:
	
	variant_name_s = ""
	variant_name_l = []
	TSSTES_s = ""
	TSSTES_l = []
	
#	print("gene_l: ", gene_l)
	#gene_l:  ['FP1', 'MYL6B', 'MYL6']	

	for i in range(len(gene_l)):
		if gene_l[i] == "TRAV14":#refの遺伝子名に"/"が含まれる
			gene_l[i] = "TRAV14/DV4"
		elif gene_l[i] == "THRA1":
			gene_l[i] = "THRA1/BTR"
		elif not gene_l[i] in isoNum_uniq_dict:
			continue
		
		if gene_l[i] == 'FP1' or gene_l[i] == 'RP1' or gene_l[i] == 'FP2' or gene_l[i] == 'RP2' or gene_l[i] == 'P1' or gene_l[i] == 'P2':
			variant_name_l.append('PRIMER')
			TSSTES_l.append('*')
			continue	
		
		if gene_l[i] == 'novel':
			variant_name_l.append('UNANNOT')
			TSSTES_l.append('*')
			continue

		else:
			if exonNum_l[i] in ref_exonNum_l[i]:
				variant_name_l.append(convert_gene_into_variant_dict[gene_l[i], exonNum_l[i]])
				TSSTES_l.append('TSS+TES')
				continue
#			print("ref_exonNum_l[i]: ", ref_exonNum_l[i])
			exonNum_l_l = exonNum_l[i].split(",")
			for j in range(len(ref_exonNum_l[i])):
#				print("j: ", j)
#				print("exonNum_l[i]: ", exonNum_l[i])
#				print("gene_l[i]: ", gene_l[i])
#				print("exonNum_l_l: ", exonNum_l_l)
#				print("ref_exonNum_l[i]: ", ref_exonNum_l[i])
#				print("ref_exonNum_l[i][j]: ", ref_exonNum_l[i][j])

				ref_exonNum_l_l = ref_exonNum_l[i][j].split(",")
#				print("ref_exonNum_l_l: ", ref_exonNum_l_l)
				eva_match = False
				if set(exonNum_l_l) <= set(ref_exonNum_l_l):
					match_index_l = []
					for k in range(len(ref_exonNum_l_l)):
						if ref_exonNum_l_l[k] in exonNum_l_l:
							match_index_l.append(k)
#					print("match_index_l: ", match_index_l, sep="\t")
					if len(match_index_l) == 1:		
						eva_match = True	
					else:
						for k in range(1, len(match_index_l)):
							if int(match_index_l[k]) == int(match_index_l[k-1]) + 1:
								eva_match = True
							else:
								eva_match = False
								break
#				print(">>> eva_match: ", eva_match, sep="\t") 
	
#				if set(exonNum_l_l) <= set(ref_exonNum_l_l):#181130変更
#				if exonNum_s in ref_exonNum_s:#181130変更
				if eva_match == True:
#					print(gene_l[i], ref_exonNum_l_l)
#					print("exonNum_s: ", exonNum_s, sep="\t")	
#					print("ref_econNum_s: ", ref_exonNum_s, sep="\t")
#					print("exonNum_l_l: ", exonNum_l_l, sep="\t")
#					print("ref_exonNum_l_l: ", ref_exonNum_l_l, sep="\t")
					if exonNum_l_l[0] == ref_exonNum_l_l[0] and exonNum_l_l[-1] != ref_exonNum_l_l[-1]:
						variant_name_s += convert_gene_into_variant_dict[gene_l[i], ref_exonNum_l[i][j]] + ","
						if variant_strand_direction_dict[convert_gene_into_variant_dict[gene_l[i], ref_exonNum_l[i][j]]] == "+":
							TSSTES_s += "TSS,"
						else:
							TSSTES_s += "TES,"
				
					elif exonNum_l_l[0] != ref_exonNum_l_l[0] and exonNum_l_l[-1] == ref_exonNum_l_l[-1]:
						variant_name_s += convert_gene_into_variant_dict[gene_l[i], ref_exonNum_l[i][j]] + ","
						if variant_strand_direction_dict[convert_gene_into_variant_dict[gene_l[i], ref_exonNum_l[i][j]]] == "+":
							TSSTES_s += "TES,"
						else:
							TSSTES_s += "TSS,"
					else:
						variant_name_s += convert_gene_into_variant_dict[gene_l[i], ref_exonNum_l[i][j]] + ","
						TSSTES_s += "None,"
				else:
##					if j == 0:
					if j == len(ref_exonNum_l[i]) -1 and variant_name_s == "":
						variant_name_s += "Novel,"
						TSSTES_s += "*,"
#						break
#					else:
#						break

			if not variant_name_s == "":
				variant_name_l.append(variant_name_s[:-1])			
				TSSTES_l.append(TSSTES_s[:-1])
				variant_name_s = ""
				TSSTES_s = ""
#	print("variant_name_l: ", variant_name_l)
#	print("TSSTES_l: ", TSSTES_l)
	
	print(line, "/".join(variant_name_l), "/".join(TSSTES_l), sep="\t")
#	print("")

