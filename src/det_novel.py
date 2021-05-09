import sys
import re

def del_terminal_index(l):
	if l[-1] == "":
		del l[-1]
	return l

f = open(sys.argv[1])
exonNum_pos_dict = {}
gene_exonNum_dict = {}
for line in f:
	line = line.replace("\n", "")
	line_l = line.split("\t")
#	print(line)
#	print(line_l[12], line_l[16], line_l[9], line_l[10], sep="\t")
	exon_num_l = line_l[16].split(",")
	startpos_l = line_l[9].split(",")
	endpos_l = line_l[10].split(",")
#	print(exon_num_l, startpos_l, endpos_l, sep="\t")
	for i in range(len(exon_num_l)):
		if not line_l[12] +":"+ exon_num_l[i] in exonNum_pos_dict:
			exonNum_pos_dict[line_l[12] +";"+ exon_num_l[i]] = [startpos_l[i], endpos_l[i]]	
#			print(exonNum_pos_dict)	
	
	if not line_l[12] in gene_exonNum_dict:
		gene_exonNum_dict[line_l[12]] = []
		gene_exonNum_dict[line_l[12]].append(line_l[16])
	else:
		gene_exonNum_dict[line_l[12]].append(line_l[16])

#for k,v in exonNum_pos_dict.items():
#	print(k,v,sep="\t")

pattern_1or2 = 0
pattern_3or4 = 0
#for k,v in gene_exonNum_dict.items():
#	print(k,v,sep="\t")
#A1BG    ['3,7,9,11,14,16,17,19', '2,7,9,11,14,16,17,19', '10,
#A1BG-AS1        ['6,9,10,17', '4,8,10,16', '2,8,12', '3,9,13'

f2 = open(sys.argv[2])
for line in f2:
	line = line.replace("\n", "")	
	line_l = line.split("\t")
#	if line_l[5] == "KNOWN":
#		continue
	if line_l[2] == "NC":
		continue
##	print(">>: ", line)
	exonNum_l = line_l[7].split(",")
	mid_exonNum_l = exonNum_l[1:-1]	
	mid_exonNum_s = ",".join(mid_exonNum_l)
#	print("exonNum_l: ", exonNum_l)
##	print("mid_exonNum_l: ", mid_exonNum_l)
	eva_exon_len_l = del_terminal_index(line_l[8].split(","))
	mid_eva_exon_len_l = eva_exon_len_l[1:-1]
#	print("eva_exon_len_l: ", eva_exon_len_l)
##	print("mid_eva_exon_len_l: ", mid_eva_exon_len_l) 

#	eva_match = False
#	for i in range(len(gene_exonNum_dict[line_l[2]])):
#		if mid_exonNum_s in gene_exonNum_dict[line_l[2]][i]:
#			print("mid_exonNum_s: ", mid_exonNum_s, sep="\t")
#			print("ref_exonNum: ", gene_exonNum_dict[line_l[2]][i], sep="\t")
#			eva_match = True
		
	eva_match = False
	eva_match_pos = ""
	if mid_exonNum_l == ['']:
		if not line_l[9] == "*":
##			print("2exon!")
			for i in range(len(gene_exonNum_dict[line_l[2]])):
				if eva_match_pos == "MATCH1MATCH2":	
					continue
				
				ref_exonNum_l = gene_exonNum_dict[line_l[2]][i].split(",")
##				print("ref_exonNum_l: ", ref_exonNum_l, sep="\t")
				
				for j in range(len(ref_exonNum_l)):
##					print(ref_exonNum_l[j], exonNum_pos_dict[line_l[2]+";"+ref_exonNum_l[j]], sep="\t")
					if exonNum_pos_dict[line_l[2]+";"+ref_exonNum_l[j]][1] == line_l[5]:
#						print("MATCH1!")
						eva_match_pos += "MATCH1"
					elif str(int(exonNum_pos_dict[line_l[2]+";"+ref_exonNum_l[j]][0])+1) == line_l[6]:
#						print("MATCH2!")
						eva_match_pos += "MATCH2"
						if eva_match_pos == "MATCH1MATCH2":
							eva_match = True
							continue
						else:
							eva_match_pos = ""
					else:
						if eva_match_pos == "MATCH1MATCH2":
							eva_match = True
							continue
						else:
							eva_match_pos = "" 
##				print("eva_match_pos: ", eva_match_pos, sep="\t")
##				print("eva_match: ", eva_match ,sep="\t")
#					continue
#			if eva_match_pos == "MATCH1MATCH2":	
#				continue
			
		else:
##			print("空")
			eva_match = True
			continue
	
	if eva_match == True:
##		print("eva_match_continue")
		continue
		
	if "short" in mid_eva_exon_len_l or "long" in mid_eva_exon_len_l or "slong" in mid_eva_exon_len_l:
		eva_match = False
	else:
		for i in range(len(gene_exonNum_dict[line_l[2]])):
			ref_exonNum_l = gene_exonNum_dict[line_l[2]][i].split(",")
##			print("ref_exonNum_l: ", ref_exonNum_l, sep="\t")
			
##			for j in range(len(ref_exonNum_l)):
##				print(ref_exonNum_l[j], exonNum_pos_dict[line_l[2]+";"+ref_exonNum_l[j]], sep="\t")
			
			if set(mid_exonNum_l) <= set(ref_exonNum_l):
				match_index_l = []
				for k in range(len(ref_exonNum_l)):
					if ref_exonNum_l[k] in mid_exonNum_l:
						match_index_l.append(k)

##				print("match_index_l: ", match_index_l, sep="\t")
				if len(match_index_l) == 1:
##					print(exonNum_pos_dict[line_l[2] +";"+ ref_exonNum_l[match_index_l[0]]])
#					print(ref_exonNum_l[::match_index_l[0]])
#					print(ref_exonNum_l[match_index_l[0]::])
					
					match_count = 0
					if match_index_l[0] == 0:
						match_count = 0
					elif len(ref_exonNum_l[::match_index_l[0]]) >= 2:	
##						print(exonNum_pos_dict[line_l[2] +";"+ ref_exonNum_l[match_index_l[0]-1]][1])
#						if int(exonNum_pos_dict[line_l[2] +";"+ ref_exonNum_l[match_index_l[0]-1]][1]) -int(sys.argv[5]) <= int(line_l[5]) <= int(exonNum_pos_dict[line_l[2] +";"+ ref_exonNum_l[match_index_l[0]-1]][1]):
						if int(exonNum_pos_dict[line_l[2] +";"+ ref_exonNum_l[match_index_l[0]-1]][1]) == int(line_l[5]):
							match_count += 1
					if len(ref_exonNum_l[match_index_l[0]::]) >= 2:
##						print(int(exonNum_pos_dict[line_l[2] +";"+ ref_exonNum_l[match_index_l[0]+1]][0])+1)
#						if int(exonNum_pos_dict[line_l[2] +";"+ ref_exonNum_l[match_index_l[0]+1]][0])+1 <= int(line_l[6]) <= int(exonNum_pos_dict[line_l[2] +";"+ ref_exonNum_l[match_index_l[0]+1]][0])+1+int(sys.argv[5]):
						if int(exonNum_pos_dict[line_l[2] +";"+ ref_exonNum_l[match_index_l[0]+1]][0])+1 == int(line_l[6]):
							match_count += 1
					
					if match_count == 2:
						eva_match = True
					else:
						eva_match = False
	
				elif len(match_index_l) == len(mid_exonNum_l):
##					print("match_index_l: ", match_index_l, sep="\t")
					for k in range(1, len(match_index_l)):
						if int(match_index_l[k]) == int(match_index_l[k-1])+1:
							eva_match = True
						else:
							eva_match = False
							break
				
					if eva_match == True:	
						match_count = 0
						if match_index_l[0] == 0:
							match_count = 0
						elif len(ref_exonNum_l[::match_index_l[0]]) >= 2:
##							print(exonNum_pos_dict[line_l[2] +";"+ ref_exonNum_l[match_index_l[0]-1]][1])
#							if int(exonNum_pos_dict[line_l[2] +";"+ ref_exonNum_l[match_index_l[0]-1]][1]) -int(sys.argv[5]) <= int(line_l[5]) <= int(exonNum_pos_dict[line_l[2] +";"+ ref_exonNum_l[match_index_l[0]-1]][1]):
							if int(exonNum_pos_dict[line_l[2] +";"+ ref_exonNum_l[match_index_l[0]-1]][1]) == int(line_l[5]):
								match_count += 1
						if len(ref_exonNum_l[match_index_l[-1]::]) >= 2:
##							print(int(exonNum_pos_dict[line_l[2] +";"+ ref_exonNum_l[match_index_l[-1]+1]][0])+1)
#							if int(exonNum_pos_dict[line_l[2] +";"+ ref_exonNum_l[match_index_l[-1]+1]][0])+1 <= int(line_l[6]) <= int(exonNum_pos_dict[line_l[2] +";"+ ref_exonNum_l[match_index_l[-1]+1]][0])+1+int(sys.argv[5]):
							if int(exonNum_pos_dict[line_l[2] +";"+ ref_exonNum_l[match_index_l[-1]+1]][0])+1 == int(line_l[6]):
								match_count += 1

						if match_count == 2:
							eva_match = True
						else:
							eva_match = False
						
##			print("eva_match: ", eva_match, sep="\t")	
			if eva_match == True:
				break
				
#	print("eva_match: ", eva_match, sep="\t")
	if not "short" in mid_eva_exon_len_l and not "long" in mid_eva_exon_len_l and not "slong" in mid_eva_exon_len_l:
		if eva_match == True:
			pattern_3or4 += 1
#			print("pattern_3or4: ", line)
			if str(sys.argv[3]) == "end":
				print(line)
		else:#中間exonの新規組み合わせ
			pattern_1or2 += 1
#			print("pattern_1or2: ", line)
			if str(sys.argv[3]) == "mid_comb":
				print(line)
	
	else:#中間exonにlong,short,slongあり
		pattern_1or2 += 1
		if str(sys.argv[3]) == "mid_sl":
			print(line)
##	print()
#		print("pattern_1or2: ", line)
#print("pattern_1or2", pattern_1or2, sep="\t")
#print("pattern_3or4", pattern_3or4, sep="\t")
