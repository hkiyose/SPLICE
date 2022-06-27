import sys
import re

def det_novel_tes(f):

	def del_terminal_index(l):
		if l[-1] == "":
			del l[-1]
		return l

	f1 = open(f)
	gene = ""
	gene_num_dict = {}
	variant_geneNum_readLen_dict = {}#variant_geneNum_readLen_dict[A2M_9116157_9067707_*,7,10,11,*_*,known,known,known,*] = [1, a2550790-5011-4638-8078-951b86eb88f3]
	junction_match_num_dict = {}
	genome_transcriptome_match_num_dict = {}
	for line in f1:
		line = line.replace("\n", "")
		line_l = line.split("\t")
		
		if line_l[2] in gene_num_dict:
			gene_num_dict[line_l[2]] += 1
		else:
			gene_num_dict[line_l[2]] = 1
			
		junction_match_count_l = line_l[9].split(",")

		if gene == "":
			gene = line_l[2]
		
		gap_l = line_l[-1].split(",")
			
		if gene == line_l[2]:
		
			k = gene + ";" + "*" +";"+ line_l[4] +";"+ line_l[5] +";"+ line_l[6] + ";" + line_l[7] + ";" + line_l[8] + ";" + ",".join(gap_l)
			
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
			k = gene + ";" + "*" +";"+ line_l[4] +";"+ line_l[5] +";"+ line_l[6] + ";" + line_l[7] + ";" + line_l[8] +";"+ ",".join(gap_l)
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

		
def det_novel_tes2(f, f_2):
	def del_terminal_index(l):
		if l[-1] == "":
			del l[-1]
		return l

	f1 = open(f)
	exonNum_pos_dict = {}
	gene_exonNum_dict = {}
	gene_exonNum_variant_dict = {}
	for line in f1:
		line = line.replace("\n", "")
		line_l = line.split("\t")
		exon_num_l = line_l[16].split(",")
		startpos_l = line_l[9].split(",")
		endpos_l = line_l[10].split(",")
		for i in range(len(exon_num_l)):
			if not line_l[12] +":"+ exon_num_l[i] in exonNum_pos_dict:
				exonNum_pos_dict[line_l[12] +";"+ exon_num_l[i]] = [startpos_l[i], endpos_l[i]]	
		if not line_l[12] in gene_exonNum_dict:
			gene_exonNum_dict[line_l[12]] = []
			gene_exonNum_dict[line_l[12]].append(line_l[16])
		else:
			gene_exonNum_dict[line_l[12]].append(line_l[16])
		
		gene_exonNum_variant_dict[line_l[12]+"/"+line_l[16]] = line_l[1]
		
	pattern_1or2 = 0
	pattern_3or4 = 0

	f2 = open(f_2)
	for line in f2:
		line = line.replace("\n", "")	
		line_l = line.split("\t")
		if line_l[2] == "NC":#181121確認
			continue
	#	print(">>: ", line)
		exonNum_l = line_l[7].split(",")
		mid_exonNum_l = exonNum_l[1:-1]	
		mid_exonNum_s = ",".join(mid_exonNum_l)
	#	print("mid_exonNum_l: ", mid_exonNum_l)
		eva_exon_len_l = del_terminal_index(line_l[8].split(","))
		mid_eva_exon_len_l = eva_exon_len_l[1:-1]
	#	print("mid_eva_exon_len_l: ", mid_eva_exon_len_l) 

		eva_match = False
		eva_match_pos = ""
		match_variant_l = []
		if mid_exonNum_l == ['']:
			if not line_l[9] == "*":
	#			print("2exon!")
	#			print("len(gene_exonNum_dict[line_l[2]]): ", len(gene_exonNum_dict[line_l[2]]), sep="\t")
				for i in range(len(gene_exonNum_dict[line_l[2]])):
					if eva_match_pos == "MATCH1MATCH2":	
						continue
					
					ref_exonNum_l = gene_exonNum_dict[line_l[2]][i].split(",")
	#				print("ref_exonNum_l: ", ref_exonNum_l, "gene_exonNum_variant_dict: " , gene_exonNum_variant_dict[line_l[2]+"/"+",".join(ref_exonNum_l)], sep="\t")
					
					for j in range(len(ref_exonNum_l)):
	#					print(ref_exonNum_l[j], exonNum_pos_dict[line_l[2]+";"+ref_exonNum_l[j]], sep="\t")
						if exonNum_pos_dict[line_l[2]+";"+ref_exonNum_l[j]][1] == line_l[5]:
							eva_match_pos += "MATCH1"
						elif str(int(exonNum_pos_dict[line_l[2]+";"+ref_exonNum_l[j]][0])+1) == line_l[6]:
							eva_match_pos += "MATCH2"
							if eva_match_pos == "MATCH1MATCH2":
								eva_match = True
	#							print("MATCH: ", gene_exonNum_variant_dict[line_l[2]+"/"+",".join(ref_exonNum_l)])
								match_variant_l.append(gene_exonNum_variant_dict[line_l[2]+"/"+",".join(ref_exonNum_l)])
								eva_match_pos = ""
	#							continue
							else:
								eva_match_pos = ""
						else:
							if eva_match_pos == "MATCH1MATCH2":
								eva_match = True
	#							print("MATCH: ", gene_exonNum_variant_dict[line_l[2]+"/"+",".join(ref_exonNum_l)])
								match_variant_l.append(gene_exonNum_variant_dict[line_l[2]+"/"+",".join(ref_exonNum_l)])
								eva_match_pos = ""
	#							continue
							else:
								eva_match_pos = "" 
	#				print("eva_match: ", eva_match ,sep="\t")
	#				print("match_variant_l: ", match_variant_l, sep="\t")
	#				print()
				
			else:
	#			print("空")
				print(line)
				eva_match = True
	#			print()
				continue
		
		if len(match_variant_l) >= 1:
			print(line_l[0], line_l[1], line_l[2], ",".join(match_variant_l), line_l[4], line_l[5], line_l[6], line_l[7], line_l[8], line_l[9], line_l[10], line_l[11], line_l[12], line_l[13], sep="\t")
	#		print("match_variant_l: ", match_variant_l, sep="\t")
	#		print()
			continue
			
		if "short" in mid_eva_exon_len_l or "long" in mid_eva_exon_len_l or "slong" in mid_eva_exon_len_l:
			eva_match = False
			continue
		else:
	#		print("len(gene_exonNum_dict[line_l[2]]): ", len(gene_exonNum_dict[line_l[2]]), sep="\t")
			for i in range(len(gene_exonNum_dict[line_l[2]])):
				ref_exonNum_l = gene_exonNum_dict[line_l[2]][i].split(",")
	#			print("ref_exonNum_l: ", ref_exonNum_l, "gene_exonNum_variant_dict: " ,gene_exonNum_variant_dict[line_l[2]+"/"+",".join(ref_exonNum_l)], sep="\t")
				
	#			for j in range(len(ref_exonNum_l)):
	#				print(ref_exonNum_l[j], exonNum_pos_dict[line_l[2]+";"+ref_exonNum_l[j]], sep="\t")
				
				if set(mid_exonNum_l) <= set(ref_exonNum_l):
					match_index_l = []
					for k in range(len(ref_exonNum_l)):
						if ref_exonNum_l[k] in mid_exonNum_l:
							match_index_l.append(k)

	#				print("match_index_l: ", match_index_l, sep="\t")
					if len(match_index_l) == 1:
	#					print(exonNum_pos_dict[line_l[2] +";"+ ref_exonNum_l[match_index_l[0]]])
						
						match_count = 0
						if match_index_l[0] == 0:
							match_count = 0
						elif len(ref_exonNum_l[::match_index_l[0]]) >= 2:	
	#						print(exonNum_pos_dict[line_l[2] +";"+ ref_exonNum_l[match_index_l[0]-1]][1])
							if int(exonNum_pos_dict[line_l[2] +";"+ ref_exonNum_l[match_index_l[0]-1]][1]) == int(line_l[5]):
								match_count += 1
						if len(ref_exonNum_l[match_index_l[0]::]) >= 2:
	#						print(int(exonNum_pos_dict[line_l[2] +";"+ ref_exonNum_l[match_index_l[0]+1]][0])+1)
							if int(exonNum_pos_dict[line_l[2] +";"+ ref_exonNum_l[match_index_l[0]+1]][0])+1 == int(line_l[6]):
								match_count += 1
						
						if match_count == 2:
							eva_match = True	
	#						print("MATCH: ", gene_exonNum_variant_dict[line_l[2]+"/"+",".join(ref_exonNum_l)])
							match_variant_l.append(gene_exonNum_variant_dict[line_l[2]+"/"+",".join(ref_exonNum_l)])
						else:
							eva_match = False
		
					elif len(match_index_l) == len(mid_exonNum_l):
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
	#							print(exonNum_pos_dict[line_l[2] +";"+ ref_exonNum_l[match_index_l[0]-1]][1])
								if int(exonNum_pos_dict[line_l[2] +";"+ ref_exonNum_l[match_index_l[0]-1]][1]) == int(line_l[5]):
									match_count += 1
							if len(ref_exonNum_l[match_index_l[-1]::]) >= 2:
	#							print(int(exonNum_pos_dict[line_l[2] +";"+ ref_exonNum_l[match_index_l[-1]+1]][0])+1)
								if int(exonNum_pos_dict[line_l[2] +";"+ ref_exonNum_l[match_index_l[-1]+1]][0])+1 == int(line_l[6]):
									match_count += 1

							if match_count == 2:
								eva_match = True
	#							print("MATCH: ", gene_exonNum_variant_dict[line_l[2]+"/"+",".join(ref_exonNum_l)])
								match_variant_l.append(gene_exonNum_variant_dict[line_l[2]+"/"+",".join(ref_exonNum_l)])
							else:
								eva_match = False
							
	#			print("eva_match: ", eva_match, sep="\t")
	#			print("match_variant_l: ", match_variant_l, sep="\t")
	#			print()	
				
		if len(match_variant_l) >= 1:
			print(line_l[0], line_l[1], line_l[2], ",".join(match_variant_l), line_l[4], line_l[5], line_l[6], line_l[7], line_l[8], line_l[9], line_l[10], line_l[11], line_l[12], line_l[13], sep="\t")
	#		print("match_variant_l: ", match_variant_l, sep="\t")


def det_novel_tes3(f, f_2):
	f1 = open(f)
	read_info_dict = {}
	for line in f1:
		line = line.replace("\n","")
		line_l = line.split("\t")
#		print(line)
		read_info_dict[line_l[0]] = ";".join(line_l)

	f2 = open(f_2)
	for line in f2:
		line = line.replace("\n","")
		line_l = line.split("\t")
#		print(line)
		read_l = line_l[-1].split(",")
		for i in range(len(read_l)):
			read_info_l = read_info_dict[read_l[i]].split(";")
#			print(read_l[i], read_info_l, sep="\t")
		
			if line_l[3] == "*":
				print("1", "0.00", line_l[2], read_info_l[14], line_l[5], line_l[6], read_info_l[12], read_info_l[13], line_l[12], read_l[i], sep="\t")
			else:
				print("1", "0.00", line_l[2], line_l[3], line_l[5], line_l[6], read_info_l[12], read_info_l[13], line_l[12], read_l[i], sep="\t")
	
def det_novel_tes4(f):

	def del_terminal_index(l):
		if l[-1] == "":
			del l[-1]
		return l

	f1 = open(f)
	gene = ""
	variant_geneNum_readLen_dict = {}#variant_geneNum_readLen_dict[A2M_9116157_9067707_*,7,10,11,*_*,known,known,known,*] = [1, a2550790-5011-4638-8078-951b86eb88f3]
	for line in f1:
		line = line.replace("\n", "")
		line_l = line.split("\t")
		
#		print(line)
		if gene == "":
			gene = line_l[2]
		if gene == line_l[2]:
			k = gene + ";" + line_l[3] +";"+ line_l[4] +";"+ line_l[5] +";"+ line_l[6] + ";" + line_l[7] +";"+ line_l[8]  
			if k in variant_geneNum_readLen_dict:
				variant_geneNum_readLen_dict[k][0] += 1
				variant_geneNum_readLen_dict[k][1] += ","+line_l[-1]
			else:
				variant_geneNum_readLen_dict[k] = [1, line_l[-1]]
		else:
			for k, v in variant_geneNum_readLen_dict.items():
				k_l = k.split(";")
				
				print(v[0], k_l[0], k_l[1], k_l[2], k_l[3], k_l[4], k_l[5], k_l[6], v[1], sep="\t")
			
			variant_geneNum_readLen_dict = {}
			
			gene = line_l[2]
			k = gene + ";" + line_l[3] +";"+ line_l[4] +";"+ line_l[5] +";"+ line_l[6] + ";" + line_l[7] +";"+ line_l[8]
			variant_geneNum_readLen_dict[k] = [1, line_l[-1]]
	
	for k, v in variant_geneNum_readLen_dict.items():
		k_l = k.split(";")
		print(v[0], k_l[0], k_l[1], k_l[2], k_l[3], k_l[4], k_l[5], k_l[6], v[1], sep="\t")

def det_novel_tes5(f, f_2):
	def del_terminal_index(l):
		if l[-1] == "":
			del l[-1]
		return l

	f1 = open(f)
	gene_exonNum_dict = {}
	gene_exonNum_variant_dict = {}
	for line in f1:
		line = line.replace("\n", "")
		line_l = line.split("\t")
		if not line_l[12] in gene_exonNum_dict:
			gene_exonNum_dict[line_l[12]] = []
			gene_exonNum_dict[line_l[12]].append(line_l[16])
		else:
			gene_exonNum_dict[line_l[12]].append(line_l[16])
		
		gene_exonNum_variant_dict[line_l[12] +"/"+ line_l[16]] = line_l[1] 

	f2 = open(f_2)
	for line in f2:
		line = line.replace("\n", "")	
		line_l = line.split("\t")
		eva_exon_len_l = del_terminal_index(line_l[6].split(","))
		exonNum_l = line_l[5].split(",")
		match_variant = ""
	
		for i in range(len(gene_exonNum_dict[line_l[1]])):
			ref_exonNum_l = gene_exonNum_dict[line_l[1]][i].split(",")
			if line_l[5] == ",".join(ref_exonNum_l):
				if not "short" in line_l[6] and not "long" in line_l[6]:
					match_variant = gene_exonNum_variant_dict[line_l[1] +"/"+ ",".join(ref_exonNum_l)]
					break
				
		if match_variant == "":
			print("NOVEL: ", line, sep="\t")
		else:	
			print("KNOWN: ", line_l[0], line_l[1], match_variant, line_l[3], line_l[4], line_l[5], line_l[6], line_l[7], line_l[8], sep="\t")

def det_novel_tes6(f, f_2, f_3, int4, int5, f_4):
	base_convert_dict = {"A":"T", "T":"A", "G":"C", "C":"G"}
	def get_rev_seq(seq):
		seq2 = ""
		for i in range(len(seq)):
			seq2 += base_convert_dict[seq[i]]
		seq2 = seq2[::-1]
		return seq2
		

	f1 = open(f)#refseq_exonNum
	gene_strand_dict = {}
	for line in f1:
		line = line.replace("\n", "")
		line_l = line.split("\t")
		if not line_l[12] in gene_strand_dict:
			gene_strand_dict[line_l[12]] = line_l[3]

	f2 = open(f_2)#
	read_3prime_seq_dict = {}#read_3prime_seq_dict[read_id] = start_pos, end_posi, strand...Invert the array when "-".
	for line in f2:
		line = line.replace("\n", "")
		line_l = line.split("\t")
		if "/" in line_l[11] or "novel" in line_l[11]:
			continue
		start_pos_l = line_l[4].split(",")
		end_pos_l = line_l[5].split(",")	
		
		if gene_strand_dict[line_l[11]] == "+":
			v = [str(int(end_pos_l[-1])+1), line_l[1], line_l[8]]
			read_3prime_seq_dict[line_l[0]] = v
			
		else:
			v = ["0", str(int(start_pos_l[0])-1), line_l[8]]
			read_3prime_seq_dict[line_l[0]] = v
		
	f3 = open(f_3)#fastq
	count = 0
	pA_count_dict = {}
	pA_read_dict = {}
	for line in f3:
		line = line.replace("\n", "")
		count += 1
		if count % 4 == 1:
			line_l = line.split()
			read_id = line_l[0][1:]
	#		print("read_id: ", read_id, sep="\t")
		
		elif count % 4 == 2:
			if read_id in read_3prime_seq_dict:
	#			print("read_3prime_seq_dict[read_id]: ", read_3prime_seq_dict[read_id], sep="\t")
				sc_3prime_seq = ""
				if read_3prime_seq_dict[read_id][2] == "+":
	#				print(line[int(read_3prime_seq_dict[read_id][0]):int(read_3prime_seq_dict[read_id][1])])
					sc_3prime_seq = line[int(read_3prime_seq_dict[read_id][0]):int(read_3prime_seq_dict[read_id][1])]
				else:
					rev_seq = get_rev_seq(line)
	#				print(rev_seq[int(read_3prime_seq_dict[read_id][0]):int(read_3prime_seq_dict[read_id][1])])
					sc_3prime_seq = rev_seq[int(read_3prime_seq_dict[read_id][0]):int(read_3prime_seq_dict[read_id][1])]

				pA_eva = False
				for i in range(len(sc_3prime_seq)-int(int4)+1):
					seg_seq = sc_3prime_seq[i:i+int(int4)]	
					if read_3prime_seq_dict[read_id][0] == "0":
						if seg_seq.count("T") >= int(int(int(int4) * float(int5))):
	#						print("polyT: ", seg_seq.count("T"))	
							pA_eva = True
					else:		
						if seg_seq.count("A") >= int(int(int(int4) * float(int5))):
	#						print("polyA: ", seg_seq.count("A"))
							pA_eva = True
				if pA_eva == True:
					pA_read_dict[read_id] = 0

	f4 = open(f_4)
	for line in f4:
		line = line.replace("\n","")
		line_l = line.split("\t")
		if "KNOWN" in line_l[0]:
			continue	
		read_id_l = line_l[-1].split(",")
		read_id_l2 = []
		for i in range(len(read_id_l)):
			if read_id_l[i] in pA_read_dict:
				read_id_l2.append(read_id_l[i])
		if len(read_id_l2) >= 1:
			print(line_l[0], len(read_id_l2), line_l[2], line_l[3], line_l[4], line_l[5], line_l[6], line_l[7], line_l[8], ",".join(read_id_l2), sep="\t")

def det_novel_tes7(f, f_2, f_3):
	def del_terminal_index(l):
		if l[-1] == "":
			del l[-1]
		return l

	f1 = open(f)
	gene_exonNum_dict = {}
	gene_strand_dict = {}
	for line in f1:
		line = line.replace("\n", "")
		line_l = line.split("\t")
		if not line_l[12] in gene_exonNum_dict:
			gene_exonNum_dict[line_l[12]] = []
			gene_exonNum_dict[line_l[12]].append(line_l[16])
			gene_strand_dict[line_l[12]] = line_l[3]
		else:
			gene_exonNum_dict[line_l[12]].append(line_l[16])

	f3 = open(f_3)
	annot_num_dict = {}
	for line in f3:
		line = line.replace("\n", "")
		line_l = line.split("\t")
		annot_l = del_terminal_index(line_l[11].split("/"))#gene
		for i in range(len(annot_l)):
			if annot_l[i] in annot_num_dict:
				annot_num_dict[annot_l[i]] += 1
			else:
				annot_num_dict[annot_l[i]] = 1

	f2 = open(f_2)
	for line in f2:
		line = line.replace("\n", "")	
		line_l = line.split("\t")
		exonNum_l = line_l[6].split(",")
		eva_exon_len_l = del_terminal_index(line_l[7].split(","))
	#	print("exonNum_l: ", exonNum_l, sep="\t")
	#	print("eva_exon_len_l: ", eva_exon_len_l, sep="\t")

		if len(exonNum_l) == 1:#1exonのみ
			continue
		
		if gene_strand_dict[line_l[2]] == "+":#Remove the exon on the 5' side. If there is still no correspondence, it can be judged as a new 3' side (polyA side).
			del exonNum_l[0]
			del eva_exon_len_l[0]
		else:
			del exonNum_l[-1]
			del eva_exon_len_l[-1]
		
	#	print("masked_exonNum_l: ", exonNum_l, sep="\t")
	#	print("masked_eva_exon_len_l: ", eva_exon_len_l, sep="\t")
		
		eva_match = False	
		for i in range(len(gene_exonNum_dict[line_l[2]])):
			ref_exonNum_l = gene_exonNum_dict[line_l[2]][i].split(",")
	#		print("ref_exonNum_l: ", ref_exonNum_l, sep="\t")
			if set(exonNum_l) <= set(ref_exonNum_l):
				match_index_l = []
	#			print(">>>ref_exonNum_l: ", ref_exonNum_l, sep="\t")
				eva_match = True
		
		if "short" in eva_exon_len_l or "long" in eva_exon_len_l or "slong" in eva_exon_len_l:
			eva_match = False

	#	print("eva_match: ", eva_match, sep="\t")	
			
		if eva_match == False:
			print(line_l[1], round(int(line_l[1])/annot_num_dict[line_l[2]],3), line_l[2], line_l[3], line_l[4], line_l[5], line_l[6], line_l[7], line_l[8], line_l[9], sep="\t")	
			
		
