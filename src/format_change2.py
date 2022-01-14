import sys
import re

def del_term_index(l):
	if l[-1] == "":
		del l[-1]
	return l

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


gene_pos_dict = {}#gene_pos_dict[gene/exon_num] = [start_pos, end_pos]
gene_strand_dict = {}
gene_splicing_variant_dict = {}#gene_splicing_variant_dict[gene] = [splicing_variant1, splicing_variant2, , , ]  
splicing_variant_UTR_dict = {}#splicing_variant_UTR_dict[splicing_variant] = [UTR5_start, UTR5_end, UTR3_start, UTR3_end]
splicing_variant_len_dict = {}
splicing_variant_exon_count_dict = {}
f3 = open(sys.argv[3])
for line in f3:
	line = line.replace("\n", "")
	line_l = line.split("\t")
#	print(line)
	start_pos_l = del_term_index(line_l[9].split(","))
	end_pos_l = del_term_index(line_l[10].split(","))
	exon_num_l = line_l[-1].split(",")
#	print("start_pos_l: ", start_pos_l, sep="\t")
#	print("end_pos_l: ", end_pos_l, sep="\t")
#	print("exon_num_l: ", exon_num_l, sep="\t")
	transcript_len = 0
	for i in range(len(start_pos_l)):
		transcript_len += int(end_pos_l[i]) - int(start_pos_l[i]) + 1
		k = line_l[12] +"/"+ exon_num_l[i]
		if not k in gene_pos_dict:
			gene_pos_dict[k] = [start_pos_l[i], end_pos_l[i]]
	
	splicing_variant_len_dict[line_l[1]] = transcript_len
	splicing_variant_exon_count_dict[line_l[1]] = len(start_pos_l)
	
#	print("transcript_len: ", transcript_len, sep="\t")	
	if not line_l[12] in gene_splicing_variant_dict:
		gene_splicing_variant_dict[line_l[12]] = []
		gene_splicing_variant_dict[line_l[12]].append(line_l[1])
	else:
		gene_splicing_variant_dict[line_l[12]].append(line_l[1])
#	print("gene_splicing_variant_dict: ", line_l[12], gene_splicing_variant_dict[line_l[12]], sep="\t")	
#	print(line_l[13],line_l[14])
	if line_l[13] == "none" or line_l[13] == "incmpl" or line_l[14] == "incmpl":
		UTR5_l = [0, 0]	
		UTR3_l = [0, 0]
	else:
		if line_l[3] == "+":
			UTR5_l = [int(line_l[4])+1, int(line_l[6])]
			UTR3_l = [int(line_l[7])+1, int(line_l[5])]
		else:
			UTR5_l = [int(line_l[7])+1, int(line_l[5])]
			UTR3_l = [int(line_l[4])+1, int(line_l[6])]	
#	print("UTR5_l: ", UTR5_l, sep="\t")
#	print("UTR3_l: ", UTR3_l, sep="\t")	
	splicing_variant_UTR_dict[line_l[1]] = [UTR5_l[0], UTR5_l[1], UTR3_l[0], UTR3_l[1]]
#	print("splicing_variant_UTR_dict: ", line_l[1], splicing_variant_UTR_dict[line_l[1]], sep="\t")
	if not line_l[12] in gene_strand_dict:
		gene_strand_dict[line_l[12]] = line_l[3]
#	print()

mt_dna_l = ["COX1", "COX2", "COX3", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "ATP6", "ATP8", "CYTB"]
sample_count = 0
error_rate_count = 0
f4 = open(sys.argv[4])
for line in f4:
	line = line.replace("\n","")
	line_l = line.split("\t")
#	print(line)
	if line.startswith("Log2"):
		for i in range(len(line_l)):
			if "_read/1Mread" in line_l[i]:
				sample_count = i+1
			if "_error_rate" in line_l[i] and error_rate_count == 0:
				error_rate_count = i
#		print("error_rate_count: ", error_rate_count ,sep="\t")
		continue
	info_l = line_l[1].split("/")
	exon_num_l = info_l[3].split(",")
	gap_l = info_l[6].split(",")
	transcript_l = info_l[1].split(",")	
	eva_exon_len_l = info_l[4].split(",")
	coding_type = "NON-CODING"
	if info_l[0] in protein_coding_dict:
		coding_type = "CODING"	
	if info_l[0] in gene_strand_dict:
		strand = gene_strand_dict[info_l[0]]
	
	if info_l[2] == "partially_known":
#		print(line)
#		print("strand: ", strand, sep="\t")
#		print("transcript_l: ", transcript_l, sep="\t")
		longest_splicing_variant = ["",0]
		longest_splicing_variant2 = ["",0]#Longest candidate regardless of UTRmatch
		for i in range(len(transcript_l)):
			if transcript_l[i] in splicing_variant_UTR_dict:
#				print(transcript_l[i], splicing_variant_len_dict[transcript_l[i]], splicing_variant_UTR_dict[transcript_l[i]], sep="\t")
				term_exon_l1 = [x for x in range(int(gene_pos_dict[info_l[0] +"/"+ exon_num_l[0]][0]) + int(gap_l[0]) +1, int(gene_pos_dict[info_l[0] +"/"+ exon_num_l[0]][1])+1)]
				term_exon_l2 = [x for x in range(int(gene_pos_dict[info_l[0] +"/"+ exon_num_l[-1]][0]) +1, int(gene_pos_dict[info_l[0] +"/"+ exon_num_l[-1]][1]) + int(gap_l[-1])+1)]
#				print("term_exon_l1: ", term_exon_l1[0], term_exon_l1[-1], sep="\t")
#				print("term_exon_l2: ", term_exon_l2[0], term_exon_l2[-1], sep="\t")
				term_exon_set1 = set(term_exon_l1)
				term_exon_set2 = set(term_exon_l2)

#				print("splicing_variant_UTR_dict[transcript_l[i]]: ", splicing_variant_UTR_dict[transcript_l[i]], sep="\t")
				UTR5_l = [x for x in range(int(splicing_variant_UTR_dict[transcript_l[i]][0]), int(splicing_variant_UTR_dict[transcript_l[i]][1])+1)]
				UTR3_l = [x for x in range(int(splicing_variant_UTR_dict[transcript_l[i]][2]), int(splicing_variant_UTR_dict[transcript_l[i]][3])+1)]
#				print("UTR5_l: ", UTR5_l[0], UTR5_l[-1], sep="\t")
#				print("UTR3_l: ", UTR3_l[0], UTR3_l[-1], sep="\t")
				UTR5_set = set(UTR5_l)
				UTR3_set = set(UTR3_l)
				
				if strand == "+":
					overlap1 = term_exon_set1 & UTR5_set
					overlap2 = term_exon_set2 & UTR3_set
				else:
					overlap1 = term_exon_set2 & UTR5_set
					overlap2 = term_exon_set1 & UTR3_set
				
				UTR5_overlap_num = len(overlap1)
				UTR3_overlap_num = len(overlap2)
#				print("UTR5_overlap_num: ", UTR5_overlap_num, sep="\t")
#				print("UTR3_overlap_num: ", UTR3_overlap_num, sep="\t")
				
				if UTR5_overlap_num >= 1 and UTR3_overlap_num >= 1 and splicing_variant_len_dict[transcript_l[i]] > longest_splicing_variant[1]:
					longest_splicing_variant = [transcript_l[i], splicing_variant_len_dict[transcript_l[i]]]
				if splicing_variant_len_dict[transcript_l[i]] > longest_splicing_variant2[1]:
					longest_splicing_variant2 = [transcript_l[i], splicing_variant_len_dict[transcript_l[i]]]
					
#		print("longest_splicing_variant: ", longest_splicing_variant, sep="\t")
#		print("longest_splicing_variant2: ", longest_splicing_variant2, sep="\t")
		if longest_splicing_variant[0] == "":
			print(coding_type, "KNOWN", "PARTIAL", info_l[0], longest_splicing_variant2[0], "-", "\t".join(line_l[2:sample_count]), sep="\t")
		else:
			print(coding_type, "KNOWN", "FULL", info_l[0], longest_splicing_variant[0], "-", "\t".join(line_l[2:sample_count]), sep="\t")
#		print()
	
	elif info_l[2] == "KNOWN":
#		print(line)
		if info_l[0] in mt_dna_l or "MT-" in info_l[0]:
			print(coding_type, "mtDNA", "FULL", info_l[0] , info_l[1], "-", "\t".join(line_l[2:sample_count]), sep="\t")
		else:
			print(coding_type, "KNOWN", "FULL", info_l[0], info_l[1], "-", "\t".join(line_l[2:sample_count]), sep="\t")
#		print()
		
	elif info_l[2] == "single_exon":
#		print(line)	
		if info_l[0] in mt_dna_l or "MT-" in info_l[0]:
			print(coding_type, "mtDNA", "PARTIAL", info_l[0] , info_l[1], "-", "\t".join(line_l[2:sample_count]), sep="\t")
		else:
#			print(">>>")
			longest_splicing_variant = ["",0,0]
			for i in range(len(transcript_l)):
				if transcript_l[i] in splicing_variant_UTR_dict:
#					print(transcript_l[i], splicing_variant_len_dict[transcript_l[i]], splicing_variant_UTR_dict[transcript_l[i]], sep="\t")
#					print("splicing_variant_exon_count_dict[transcript_l[i]]: ", splicing_variant_exon_count_dict[transcript_l[i]], sep="\t")
					if splicing_variant_len_dict[transcript_l[i]] > longest_splicing_variant[1]:
						longest_splicing_variant = [transcript_l[i], splicing_variant_len_dict[transcript_l[i]], splicing_variant_exon_count_dict[transcript_l[i]]]
#			print("longest_splicing_variant: ", longest_splicing_variant, sep="\t")	
			if longest_splicing_variant[2] == 1:#referenceがもともとsingle_exon
				print(coding_type, "KNOWN", "PARTIAL", info_l[0] , longest_splicing_variant[0], "-", "\t".join(line_l[2:sample_count]), sep="\t")	
			else:
				print(coding_type, "KNOWN", "SINGLE", info_l[0], longest_splicing_variant[0], "-", "\t".join(line_l[2:sample_count]), sep="\t")
#		print()
	elif info_l[2] == "3UTR":
#		print(line)
#		print("strand: ", strand, sep="\t")
		UTR3_exon_gap = 0
		if strand == "+":
			UTR3_exon_gap = gap_l[-1]	
		else:
			UTR3_exon_gap = gap_l[0]
#		print("UTR3_exon_gap: ", UTR3_exon_gap, sep="\t")
		longest_splicing_variant = ["",0]#Candidates with UTR-UTR and small 3'UTRgap (less than 30bp) will be preferred for long splicingVariant.=>KNOWN,FULL,GENE,VARIANT
		longest_splicing_variant2 = ["",0,1000000]#For candidates with UTR-UTR and large 3'UTRgap (>30bp), give priority to the candidate with smaller 3'UTRgap (if they are the same, give priority to the candidate with >longer 3'UTRgap)=>KNOWN,FULL_3UTR,GENE,VARIANT
		longest_splicing_variant3 = ["",0]#KNOWN, PARTIAL, GENE, VARIANT
		longest_splicing_variant4 = ["",0,1000000]#KNOWN, PARTIAL_3UTR, GENE, VARIANT
		longest_splicing_variant5 = ["",0]#KNOWN, PARTIAL, GENE, VARIANT#non-codingなど 
		for i in range(len(transcript_l)):
			if transcript_l[i] in splicing_variant_UTR_dict:
#				print(transcript_l[i], splicing_variant_len_dict[transcript_l[i]], splicing_variant_UTR_dict[transcript_l[i]], sep="\t")
				term_exon_l1 = [x for x in range(int(gene_pos_dict[info_l[0] +"/"+ exon_num_l[0]][0]) + int(gap_l[0]) +1, int(gene_pos_dict[info_l[0] +"/"+ exon_num_l[0]][1])+1)]
				term_exon_l2 = [x for x in range(int(gene_pos_dict[info_l[0] +"/"+ exon_num_l[-1]][0]) +1, int(gene_pos_dict[info_l[0] +"/"+ exon_num_l[-1]][1]) + int(gap_l[-1])+1)]
#				print("term_exon_l1: ", term_exon_l1[0], term_exon_l1[-1], sep="\t")
#				print("term_exon_l2: ", term_exon_l2[0], term_exon_l2[-1], sep="\t")		
				term_exon_set1 = set(term_exon_l1)
				term_exon_set2 = set(term_exon_l2)
	
				UTR5_l = [x for x in range(int(splicing_variant_UTR_dict[transcript_l[i]][0]), int(splicing_variant_UTR_dict[transcript_l[i]][1])+1)]
				UTR3_l = [x for x in range(int(splicing_variant_UTR_dict[transcript_l[i]][2]), int(splicing_variant_UTR_dict[transcript_l[i]][3])+1)]
#				print("UTR5_l: ", UTR5_l[0], UTR5_l[-1], sep="\t")
#				print("UTR3_l: ", UTR3_l[0], UTR3_l[-1], sep="\t")
				UTR5_set = set(UTR5_l)
				UTR3_set = set(UTR3_l)
                                
				if strand == "+":
					overlap1 = term_exon_set1 & UTR5_set
					overlap2 = term_exon_set2 & UTR3_set
				else:
					overlap1 = term_exon_set2 & UTR5_set
					overlap2 = term_exon_set1 & UTR3_set
                                
				UTR5_overlap_num = len(overlap1)
				UTR3_overlap_num = len(overlap2)
#				print("UTR5_overlap_num: ", UTR5_overlap_num, sep="\t")
#				print("UTR3_overlap_num: ", UTR3_overlap_num, sep="\t")			
				
				if strand == "+":
					UTR3_exon_gap = term_exon_l2[-1] - int(splicing_variant_UTR_dict[transcript_l[i]][3])
				else:
					UTR3_exon_gap = term_exon_l1[0] - int(splicing_variant_UTR_dict[transcript_l[i]][3])
#				print("UTR3_exon_gap: ", UTR3_exon_gap, sep="\t")
				
				if UTR5_overlap_num >= 1 and UTR3_overlap_num >= 1:#KNOWN,FULL
					if abs(UTR3_exon_gap) < int(sys.argv[5]):
						if splicing_variant_len_dict[transcript_l[i]] > longest_splicing_variant[1]:
							longest_splicing_variant = [transcript_l[i], splicing_variant_len_dict[transcript_l[i]]]
					else:#3UTR_novel	
						if abs(UTR3_exon_gap) == abs(longest_splicing_variant2[2]):
							if splicing_variant_len_dict[transcript_l[i]] > longest_splicing_variant2[1]:
								longest_splicing_variant2 = [transcript_l[i], splicing_variant_len_dict[transcript_l[i]], UTR3_exon_gap]
						elif abs(UTR3_exon_gap) < abs(longest_splicing_variant2[2]):
							longest_splicing_variant2 = [transcript_l[i], splicing_variant_len_dict[transcript_l[i]], UTR3_exon_gap]
				
				else:#KNOWN,PARTIAL
					if abs(UTR3_exon_gap) < int(sys.argv[5]):
						if splicing_variant_len_dict[transcript_l[i]] > longest_splicing_variant3[1]:
							longest_splicing_variant3 = [transcript_l[i], splicing_variant_len_dict[transcript_l[i]]]
					elif abs(UTR3_exon_gap) == abs(term_exon_l2[-1]) or abs(UTR3_exon_gap) == abs(term_exon_l1[0]):#non-codingなど
						if splicing_variant_UTR_dict[transcript_l[i]] == [0,0,0,0]:
							if splicing_variant_len_dict[transcript_l[i]] > longest_splicing_variant5[1]:
								longest_splicing_variant5 = [transcript_l[i], splicing_variant_len_dict[transcript_l[i]]]
					else:#3UTR_novel		
						if abs(UTR3_exon_gap) == abs(longest_splicing_variant4[2]):
							if splicing_variant_len_dict[transcript_l[i]] > longest_splicing_variant4[1]:
								longest_splicing_variant4 = [transcript_l[i], splicing_variant_len_dict[transcript_l[i]], UTR3_exon_gap]
						elif abs(UTR3_exon_gap) < abs(longest_splicing_variant4[2]):
								longest_splicing_variant4 = [transcript_l[i], splicing_variant_len_dict[transcript_l[i]], UTR3_exon_gap]
						
#					print(">>>")
				
#				print("--------------------")
#		print("longest_splicing_variant: ", longest_splicing_variant, sep="\t")
#		print("longest_splicing_variant2: ", longest_splicing_variant2, sep="\t")
#		print("longest_splicing_variant3: ", longest_splicing_variant3, sep="\t")
#		print("longest_splicing_variant4: ", longest_splicing_variant4, sep="\t")
#		print("longest_splicing_variant5: ", longest_splicing_variant5, sep="\t")
		if not longest_splicing_variant[0] == "":
			print(coding_type, "KNOWN", "FULL", info_l[0] , longest_splicing_variant[0], "-", "\t".join(line_l[2:sample_count]), sep="\t")
		elif not longest_splicing_variant2[0] == "":
			if strand == "+" and longest_splicing_variant2[2] > 0 or strand == "-" and longest_splicing_variant2[2] < 0:
				print(coding_type, "KNOWN", "FULL", info_l[0] , longest_splicing_variant2[0], "3UTR_LONG", "\t".join(line_l[2:sample_count]), sep="\t")	
			else:
				print(coding_type, "KNOWN", "FULL", info_l[0] , longest_splicing_variant2[0], "3UTR_SHORT", "\t".join(line_l[2:sample_count]), sep="\t")
		elif not longest_splicing_variant3[0] == "":
			print(coding_type, "KNOWN", "PARTIAL", info_l[0] , longest_splicing_variant3[0], "-", "\t".join(line_l[2:sample_count]), sep="\t")
		elif not longest_splicing_variant4[0] == "":
			if strand == "+" and longest_splicing_variant4[2] > 0 or strand == "-" and longest_splicing_variant4[2] < 0:
				print(coding_type, "KNOWN", "PARTIAL", info_l[0] , longest_splicing_variant4[0], "3UTR_LONG", "\t".join(line_l[2:sample_count]), sep="\t")
			else:
				print(coding_type, "KNOWN", "PARTIAL", info_l[0] , longest_splicing_variant4[0], "3UTR_SHORT", "\t".join(line_l[2:sample_count]), sep="\t")
		elif not longest_splicing_variant5[0] == "":
			print(coding_type, "KNOWN", "PARTIAL", info_l[0] , longest_splicing_variant5[0], "-", "\t".join(line_l[2:sample_count]), sep="\t")
#			print("ERROR2")	
#		else:
#			print("ERROR3")
		
#		print()
	else:
#		print(line)
#		print(info_l)
		max_gap = 0
		for i in range(len(gap_l)):
			if gap_l[i] == "*":
				continue
			if max_gap < abs(int(gap_l[i])):
				max_gap = abs(int(gap_l[i]))
#		print("max_gap: ", max_gap, sep="\t")
		
		error_rate_l = line_l[error_rate_count:]
		error_rate_l2 = [0]*(len(gap_l)-2)	
#		print(len(error_rate_l), error_rate_l, sep="\t")
		for i in range(len(error_rate_l)):
#			print(error_rate_l[i])
			error_rate_l_l = error_rate_l[i].split(",")
			for j in range(len(error_rate_l_l)):
				error_rate_l2[j] += int(error_rate_l_l[j])	
#		print("error_rate_l2: ", error_rate_l2, sep="\t")
		error_rate_l3 = []
		for i in range(len(error_rate_l2)):
			error_rate_l3.append(round(error_rate_l2[i]/sample_count-3,1))
#		print("error_rate_l3: ", error_rate_l3, sep="\t")
		
		max_error_rate = 0
		if max_gap == 0:#COMBINATION
			for i in range(len(error_rate_l3)):
				if int(error_rate_l3[i]) > max_error_rate:
					max_error_rate = int(error_rate_l3[i])
#			print("max_error_rate: ", max_error_rate, sep="\t")
			if max_error_rate < int(sys.argv[6]):	
				print(coding_type, "NOVEL", "COMBINATION", info_l[0], "-", "/".join(info_l[3:]), "\t".join(line_l[2:sample_count]), sep="\t")
#			else:		
#				print("ERROR")
		
		elif max_gap >= int(sys.argv[7]):
			if max_gap < int(sys.argv[8]):
				for i in range(len(eva_exon_len_l)*2):
					if i == 0 or i == len(eva_exon_len_l)*2-1:
						continue
#					print(str(i), eva_exon_len_l[int(i/2)], error_rate_l3[i-1], sep="\t")
					if not "s" in eva_exon_len_l and not "l" in eva_exon_len_l and not "sl" in eva_exon_len_l:
						if i == 1 or i == 2 or i == len(eva_exon_len_l)*2-2 or i == len(eva_exon_len_l)*2-3:
							if int(error_rate_l3[i-1]) > max_error_rate:
								max_error_rate = int(error_rate_l3[i-1])
					else:
						if not eva_exon_len_l[int(i/2)] == "k" and not eva_exon_len_l[int(i/2)] == "*":
#							print("check_error_rate!!!", error_rate_l3[i-2], error_rate_l3[i-1], error_rate_l3[i], sep="\t")
							if int(error_rate_l3[i-1]) > max_error_rate:
								max_error_rate = int(error_rate_l3[i-1])
							if int(error_rate_l3[i-2]) > max_error_rate:
								max_error_rate = int(error_rate_l3[i-2])
							if int(error_rate_l3[i]) > max_error_rate:
								max_error_rate = int(error_rate_l3[i])
#				print("max_error_rate: ", max_error_rate, sep="\t")
				if max_error_rate < int(sys.argv[9]):
					print(coding_type, "NOVEL", "LENGTH", info_l[0], "-", "/".join(info_l[3:]), "\t".join(line_l[2:sample_count]), sep="\t")
#				else:
#					print("ERROR")
			else:	
				print(coding_type, "NOVEL", "LENGTH", info_l[0], "-", "/".join(info_l[3:]), "\t".join(line_l[2:sample_count]), sep="\t")
#		else:
#			print("ERROR: ", line, sep="\t")	
			
#		print()
				
