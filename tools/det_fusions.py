import sys
import re
#import align2

def det_fusion(f, f_2, f_3, f_4, int5, int6):

	def SW(seq1, seq2, match, mismatch, gap):
		score = []
		pointer = []
		for i in range(len(seq1)):
			tmp = []
			tmp2 = []
			for j in range(len(seq2)):
				tmp.append(0)
				tmp2.append("")
			score.append(tmp)
			pointer.append(tmp2)
		max_score = 0
		max_i = 0
		max_j = 0
		for i in range(len(seq1)):
			for j in range(len(seq2)):
				if i - 1 >= 0 and j - 1 >= 0:
					(diagonal_score, left_score, up_score) = (score[i - 1][j - 1], score[i][j - 1] + gap, score[i - 1][j] + gap)
				else:
					(diagonal_score, left_score, up_score) = (0, 0 + gap, 0 + gap)

				if seq1[i] == seq2[j]:
					diagonal_score += match
				else:
					diagonal_score += mismatch
				if diagonal_score > left_score and diagonal_score > up_score:
					pointer[i][j] = "d"
					if diagonal_score > 0:
						score[i][j] = diagonal_score
					else:
						score[i][j] = 0
				elif left_score > diagonal_score and left_score > up_score:
					pointer[i][j] = "l"
					if left_score > 0:
						score[i][j] = left_score
					else:
						score[i][j] = 0
				else:
					pointer[i][j] = "u"
					score[i][j] = up_score
					if up_score > 0:
						score[i][j] = up_score
					else:
						score[i][j] = 0
				if max_score < score[i][j]:
					max_score = score[i][j]
					max_i = i
					max_j = j

		seq1_alignd = ""
		seq2_alignd = ""
		i = max_i
		j = max_j
		i2 = 0
		j2 = 0	
		while i >=0 and j >= 0 and score[i][j] > 0:
			if pointer[i][j] == "d":
				seq1_alignd += seq1[i]
				seq2_alignd += seq2[j]
				i2 += 1
				j2 += 1
				i = i - 1
				j = j - 1
			elif pointer[i][j] == "l":
				seq1_alignd += "-"
				seq2_alignd += seq2[j]
				j2 += 1
				i = i
				j = j - 1
            
			elif pointer[i][j] == "u":
				seq1_alignd += seq1[i]
				seq2_alignd += "-"
				i2 += 1
				i = i - 1
				j = j
        
		return seq1_alignd[::-1], seq2_alignd[::-1], max_i - i2 + 1, max_i, max_j - j2 + 1, max_j

	def del_terminal_index(l):
			if l[-1] == "":
					del l[-1]
			return l

	def get_primerSeq(query_seq):
			primer_seq = "AAGCAGTGGTATCAACGCAGAGTAC"#SMARTer II A Oligonucleotide
			rev_primer_seq = "GTACTCTGCGTTGATACCACTGCTT"

#			res1 = align2.SW(primer_seq, query_seq, 1, -1, -1)
#			res2 = align2.SW(rev_primer_seq, query_seq, 1, -1, -1)
			res1 = SW(primer_seq, query_seq, 1, -1, -1)
			res2 = SW(rev_primer_seq, query_seq, 1, -1, -1)

			consensus1 = ""
			match1 = 0

			for i in range(len(res1[0])):
					if res1[0][i] == res1[1][i]:
							match1 += 1
							consensus1 += "|"
					else:
							consensus1 += "*"

	##      print(res1[0])
	##      print(consensus1)
	##      print(res1[1])

			consensus2 = ""
			match2 = 0

			for i in range(len(res2[0])):
					if res2[0][i] == res2[1][i]:
							match2 += 1
							consensus2 += "|"
					else:
							consensus2 += "*"

	##      print(res2[0])
	##      print(consensus2)
	##      print(res2[1])

	##      print("match1: ", match1, match1/len(res1[0]), sep="\t")
	##      print("match2: ", match2, match2/len(res2[0]), sep="\t")

			if match1 >= int(int5) and match1/len(res1[0]) >= float(int6) or match2 >= int(int5) and match2/len(res2[0]) >= float(int6):
	##              print("MATCH!")
					return "match"
			else:
	##              print("unMATCH!")
					return "unmatch"

	f1 = open(f)
	gene_exonnum_pos_dict = {}#gene_exonnum_pos_dict["gene/exonnum"] = ["start_pos", "end_pos"]
	for line in f1:
		line = line.replace("\n","")
		line_l = line.split("\t")
		start_pos_l = line_l[9].split(",")
		end_pos_l = line_l[10].split(",")	
		exon_num_l = line_l[16].split(",")
	#	print(line)
		for i in range(len(exon_num_l)):
			k = line_l[12]+"/"+exon_num_l[i]
			if not k in gene_exonnum_pos_dict:
				gene_exonnum_pos_dict[k] = ["0","0"]
				gene_exonnum_pos_dict[k][0] = start_pos_l[i] 
				gene_exonnum_pos_dict[k][1] = end_pos_l[i]
	#	print(gene_exonnum_pos_dict)

	f3 = open(f_3)
	id_fusion_gene_dict = {}
	for line in f3:
		line = line.replace("\n","")
		line_l = line.split("\t")
		annot_gene_variant_l = line_l[3].split("/")
		if len(annot_gene_variant_l) > 2:
			count = 0
			fusion_gene_l = []
			for i in range(len(annot_gene_variant_l)):
				count += 1
				if count % 2 == 0:
					fusion_gene_l.append(annot_gene_variant_l[i])
			id_fusion_gene_dict[line_l[0]] = "/".join(sorted(fusion_gene_l))
		
	#for k,v in id_fusion_gene_dict.items():
	#	print(k,v,sep="\t")	

	f4 = open(f_4)#/home/hkiyose/nanopore_RNAseq/data/RK019_RNA_190125_FAK37469_180125.fastq.rmdup
	count = 0
	reaa_id_seq_dict = {}#reaa_id_seq_dict[readID] = seq
	for line in f4:
		line = line.replace("\n","")	
		line_l = line.split()
		count += 1
		if count % 4 == 1:  
	#		print(line_l[0][1:])
			read_id = line_l[0][1:]
		elif count % 4 == 2:
	#		print(line)
			reaa_id_seq_dict[read_id] = line

	#for k,v in reaa_id_seq_dict.items():
	#	print(k, v, sep="\t")	
		
	f2 = open(f_2)
	gene_read_num_dict = {}#gene_read_num_dict[gene] = total_read_num
	readID_gene_strand_dict = {}#readID_gene_strand_dict["readID/gene"] = strand
	fusion_table_dict = {}#fusion_table_dict["gene1/gene2/chr1/chr2/bp1/bp1"] = readID/readID ~
	###gap_len_dict = {}
	###for i in range(0,101):
	###	gap_len_dict[i*100000] = [0,0]
		
	for line in f2:
		line = line.replace("\n","")
		line_l = line.split("\t")
		gene_l = del_terminal_index(line_l[11].split("/"))
		strand_l = del_terminal_index(line_l[8].split("/"))
		chr_l = del_terminal_index(line_l[3].split("/"))
		start_pos_in_read_l = del_terminal_index(line_l[4].split("/"))
		end_pos_in_read_l = del_terminal_index(line_l[5].split("/"))
		start_pos_l = del_terminal_index(line_l[6].split("/"))
		end_pos_l = del_terminal_index(line_l[7].split("/"))
		exon_num_l = del_terminal_index(line_l[12].split("/"))
	#	print(line)
		for i in range(len(gene_l)):
			if gene_l[i] in gene_read_num_dict:
				gene_read_num_dict[gene_l[i]] += 1
			else:
				gene_read_num_dict[gene_l[i]] = 1
	#	print(gene_read_num_dict)
		if len(strand_l) == 1:
			for i in range(len(gene_l)):
				readID_gene_strand_dict[line_l[0]+"/"+gene_l[i]] = strand_l[0]
		elif len(strand_l) == len(gene_l):
			for i in range(len(gene_l)):
				readID_gene_strand_dict[line_l[0]+"/"+gene_l[i]] = strand_l[i]
	#	print(readID_gene_strand_dict)
			
		#<><><><><> Get the position (genome) across the break point <><><><><>
		if len(gene_l) == 2 and not "novel" in gene_l and len(start_pos_l) == 2:#divided by minimap2
	##		print(">>>", line, sep="\t")
			start_pos_l_gene1 = start_pos_l[0].split(",")
			end_pos_l_gene1 = end_pos_l[0].split(",")
			if readID_gene_strand_dict[line_l[0]+"/"+gene_l[0]] == "+":
				gene1_breakpoint = end_pos_l_gene1[-1]
			else:
				gene1_breakpoint = start_pos_l_gene1[0]
			
			start_pos_l_gene2 = start_pos_l[1].split(",")
			end_pos_l_gene2 = end_pos_l[1].split(",")
			if readID_gene_strand_dict[line_l[0]+"/"+gene_l[1]] == "+":
				gene2_breakpoint = start_pos_l_gene2[0]
			else:
				gene2_breakpoint = end_pos_l_gene2[-1]

	##		print("gene1_breakpoint: ", gene1_breakpoint, sep="\t")
	##		print("gene2_breakpoint: ", gene2_breakpoint, sep="\t")
				
		elif len(gene_l) == 2 and not "novel" in gene_l and len(start_pos_l) == 1:#conjoined geneなど
	##		print(">>>2", line, sep="\t")
	##		print("exon_num_l: ", exon_num_l, sep="\t")
			exon_num_l_gene1 = del_terminal_index(exon_num_l[0].split(","))
	#		print("len(exon_num_l_gene1): ", len(exon_num_l_gene1), sep="\t")	
			start_pos_l_l = start_pos_l[0].split(",")
			end_pos_l_l = end_pos_l[0].split(",")
			gene1_breakpoint = end_pos_l_l[len(exon_num_l_gene1)-1]
			gene2_breakpoint = start_pos_l_l[len(exon_num_l_gene1)]
		
	##		print("gene1_breakpoint: ", gene1_breakpoint, sep="\t")
	##		print("gene2_breakpoint: ", gene2_breakpoint, sep="\t")

		else:
			continue

		#<><><><><> Get the position(read) across the break point <><><><><>
		if len(gene_l) == 2 and not "novel" in gene_l and len(start_pos_l) == 2:#divided by minimap2
			end_pos_in_read_l_gene1 = end_pos_in_read_l[0].split(",")
			start_pos_in_read_l_gene2 = start_pos_in_read_l[1].split(",")
			
			gene1_breakpoint_in_read = end_pos_in_read_l_gene1[-1]
			gene2_breakpoint_in_read = start_pos_in_read_l_gene2[0]

	##		print("gene1_breakpoint_in_read: ", gene1_breakpoint_in_read, sep="\t") 
	##		print("gene2_breakpoint_in_read: ", gene2_breakpoint_in_read, sep="\t")
			
		elif len(gene_l) == 2 and not "novel" in gene_l and len(start_pos_l) == 1:
			start_pos_in_read_l_l = start_pos_in_read_l[0].split(",") 
			end_pos_in_read_l_l = end_pos_in_read_l[0].split(",")
			
			gene1_breakpoint_in_read = end_pos_in_read_l_l[len(exon_num_l_gene1)-1]
			gene2_breakpoint_in_read = start_pos_in_read_l_l[len(exon_num_l_gene1)]
	##		print("gene1_breakpoint_in_read: ", gene1_breakpoint_in_read, sep="\t")
	##		print("gene2_breakpoint_in_read: ", gene2_breakpoint_in_read, sep="\t")
		
		eva_artificial_fusion  = False
		
		if int(gene1_breakpoint_in_read) < int(gene2_breakpoint_in_read):
	#		print(reaa_id_seq_dict[line_l[0]][int(gene1_breakpoint_in_read)+1:int(gene2_breakpoint_in_read)])
			if len(reaa_id_seq_dict[line_l[0]][int(gene1_breakpoint_in_read)+1:int(gene2_breakpoint_in_read)]) >= 10:
				eva_match = get_primerSeq(reaa_id_seq_dict[line_l[0]][int(gene1_breakpoint_in_read)+1:int(gene2_breakpoint_in_read)])
				if eva_match == "match":
					eva_artificial_fusion = True

		if eva_artificial_fusion == True:
			continue
		
	##	print("gene_1", "gene_2", "chr_1", "chr_2", "gene1_bp", "gene2_bp", sep="\t")	
	##	if len(start_pos_l) == 1:
	##		print(gene_l[0], gene_l[1], chr_l[0], chr_l[0], gene1_breakpoint, gene2_breakpoint, sep="\t")
	##	else:		
	##		print(gene_l[0], gene_l[1], chr_l[0], chr_l[1], gene1_breakpoint, gene2_breakpoint, sep="\t")

		if line_l[0] in id_fusion_gene_dict:
			if gene_l[0] == gene_l[1]:
				continue
	##		print("id_fusion_gene_dict[line_l[0]]", id_fusion_gene_dict[line_l[0]], sep="\t")
			gene_l_sort = sorted(gene_l)
	##		print("gene_l_sort: ", gene_l_sort, sep="\t")
		
			if id_fusion_gene_dict[line_l[0]] == "/".join(gene_l_sort):
	##			print("MATCH!")
				gene_dict = {}#gene_dict[gene] = [chr, gene_bp]   
				if len(start_pos_l) == 1:
					gene_dict[gene_l[0]] = [chr_l[0], gene1_breakpoint]
					gene_dict[gene_l[1]] = [chr_l[0], gene2_breakpoint]
				else:
					gene_dict[gene_l[0]] = [chr_l[0], gene1_breakpoint]
					gene_dict[gene_l[1]] = [chr_l[1], gene2_breakpoint]
	##			print("gene_dict[gene_l[0]]: ", gene_dict[gene_l[0]], sep="\t")
	##			print("gene_dict[gene_l[1]]: ", gene_dict[gene_l[1]], sep="\t")
				
				k = gene_l_sort[0]+"/"+gene_l_sort[1]+"/"+gene_dict[gene_l_sort[0]][0]+"/"+gene_dict[gene_l_sort[1]][0]+"/"+gene_dict[gene_l_sort[0]][1]+"/"+gene_dict[gene_l_sort[1]][1]
				if k in fusion_table_dict:
					fusion_table_dict[k].append(line_l[0])
				else:
					fusion_table_dict[k] = []
					fusion_table_dict[k].append(line_l[0])
	##			print("k", k, sep="\t")
	##	print("")

	###for k,v in sorted(gap_len_dict.items(), key=lambda x:x[0]):
	###	print(k, v[0], v[1], sep="\t")

	for k,v in fusion_table_dict.items():
	#	print(k, v, sep="\t")
		k_l = k.split("/")
		gene1_read_num = gene_read_num_dict[k_l[0]]
		gene2_read_num = gene_read_num_dict[k_l[1]]
	#	print("gene1_read_num: ", gene1_read_num, sep="\t")
	#	print("gene2_read_num: ", gene2_read_num, sep="\t")
		
		print(len(v), round(len(v)/((gene1_read_num+gene2_read_num)/2)*100, 3), k_l[0]+"/"+k_l[1], k_l[2]+"/"+k_l[3], k_l[4]+"/"+k_l[5], ",".join(v), sep="\t")
		
def det_fusion2(f, int2):

	f1 = open(f)
	for line in f1:
		line = line.replace("\n","")
		line_l = line.split("\t")
#		print(line)
	
		chr_l = line_l[3].split("/")
		breakpoint_l = line_l[4].split("/")
		breakpoint_gap = abs(int(breakpoint_l[0]) - int(breakpoint_l[1]))	
#		print(line)
#		print("breakpoint_gap: ", breakpoint_gap, sep="\t")
	
		if chr_l[0] == chr_l[1]:
			if breakpoint_gap <= int(int2):
#				print(">>>close")
#				print()
				continue
		print(line)

def det_fusion3(f, int2, int3, int4):

	merged_gene_bp_pos_dict = {}#merged_gene_bp_pos_dict[gene_A/gene_B] = [gene_A_bp_start_pos1, gene_A_bp_start_pos2, gene_B_bp_start_pos1, gene_B_bp_start_pos2]
	f1 = open(f)
	for line in f1:
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
				if int(merged_gene_bp_pos_dict[k][i][0]) - int(int2) <= int(bp_pos_l[0]) <= int(merged_gene_bp_pos_dict[k][i][1]) + int(int2) and int(merged_gene_bp_pos_dict[k][i][2]) - int(int2) <= int(bp_pos_l[1]) <= int(merged_gene_bp_pos_dict[k][i][3]) + int(int2) and line_l[3] == merged_gene_bp_pos_dict[k][i][7]:
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
			if v[i][4] >= int(int3) and v[i][5] >= float(int4):
				print(v[i][4], round(v[i][5],3), k, v[i][7], v[i][0]+"-"+v[i][1]+"/"+v[i][2]+"-"+v[i][3], v[i][8], sep="\t")
