import sys
import re

def merge_transcript(f, f_2):

	def del_terminal_index(l):
		if l[-1] == "":
			del l[-1]
		return l

	f1 = open(f)
	readID_merge_pos_dict = {}
	readID_junction_match_num_dict = {}
	readID_genome_transcriptome_match_num_dict = {}
	for line in f1:
		line = line.replace("\n", "")
		line_l = line.split("\t")
		start_pos_l = line_l[6].split(",")
		end_pos_l = line_l[7].split(",")
#		print(line)
		readID_merge_pos_dict[line_l[0]] = [end_pos_l[0], start_pos_l[-1]]
		readID_junction_match_num_dict[line_l[0]] = line_l[16]
		readID_genome_transcriptome_match_num_dict[line_l[0]] = [line_l[17],line_l[18]]
#	        print(readID_merge_pos_dict[line_l[0]])
	
	f2 = open(f_2)
	for line in f2:
		line = line.replace("\n", "")
		line_l = line.split("\t")
		exon_num_l = line_l[5].split(",")	
		eva_exonlen_l = del_terminal_index(line_l[6].split(","))
		readID_l = line_l[9].split(",")
##		print(line)
	
		for i in range(len(readID_l)):
#			print("1", "0.0", line_l[2], line_l[3], line_l[4], gene_exonNum_pos_dict[line_l[2] +"," + exon_num_l[0]][1], gene_exonNum_pos_dict[line_l[2] +","+ exon_num_l[-1]][0], "*,"+",".join(exon_num_l[1:-1])+",*", "*,"+",".join(eva_exonlen_l[1:-1])+",*", line_l[9], sep="\t" )
			print("1", "0.0", line_l[2], line_l[3], line_l[4], readID_merge_pos_dict[readID_l[i]][0], readID_merge_pos_dict[readID_l[i]][1], "*,"+",".join(exon_num_l[1:-1])+",*", "*,"+",".join(eva_exonlen_l[1:-1])+",*", readID_junction_match_num_dict[readID_l[i]], readID_genome_transcriptome_match_num_dict[readID_l[i]][0], readID_genome_transcriptome_match_num_dict[readID_l[i]][1], readID_l[i], sep="\t" )
	
def merge_transcript2(f, f_2, f_3):

	def del_term_index(l):
		if l[-1] == "":
			del l[-1]
		return l

	gene_exon_num_pos_dict = {}
	f1 = open(f)
	for line in f1:
		line = line.replace("\n","")
		line_l = line.split("\t")
		exon_num_l = line_l[-1].split(",")
		start_pos_l = del_term_index(line_l[9].split(","))
		end_pos_l = del_term_index(line_l[10].split(","))	
	
		for i in range(len(exon_num_l)):
			k = line_l[12] +"/"+ exon_num_l[i]
			if not k in gene_exon_num_pos_dict:
				gene_exon_num_pos_dict[k] = (start_pos_l[i], end_pos_l[i])

#gene_exon_num_pos_dict[A1BG/3] = ('58346849', '58347029')

	read_dict = {}
	f2 = open(f_2)
	for line in f2:
		line = line.replace("\n","")
		line_l = line.split("\t")
		read_dict[line_l[-1]] = 0
	
	read_gap_dict = {}
	f3 = open(f_3)
	for line in f3:
		line = line.replace("\n","")
		line_l = line.split("\t")
		if not line_l[0] in read_dict:
			continue	
	
#		print(line_l[11], line_l[12], line_l[13], sep="\t")
#		print(line_l[6], line_l[7], sep="\t")

		exon_num_l = line_l[12].split(",")
		eva_exon_len_l = line_l[13].split(",")
		start_pos_l = line_l[6].split(",")
		end_pos_l = line_l[7].split(",")
#		print("exon_num_l: ", exon_num_l, "eva_exon_len_l: ", eva_exon_len_l, sep="\t")
#		print("start_pos_l: ", start_pos_l, "end_pos_l: ", end_pos_l, sep="\t")
	
		gap_l = []
		for i in range(len(exon_num_l)):
#			if i == 0 or i ==  len(exon_num_l)-1:
#				continue
			k = line_l[11] +"/"+ exon_num_l[i]
			
			ref_start_pos = int(gene_exon_num_pos_dict[k][0]) + 1
			ref_end_pos = int(gene_exon_num_pos_dict[k][1])
			read_start_pos = int(start_pos_l[i])
			read_end_pos = int(end_pos_l[i])
#			print("ref_start_pos: ", ref_start_pos, "ref_end_pos: ", ref_end_pos, "read_start_pos: ", read_start_pos, "read_end_pos: ", read_end_pos, sep="\t")
			
			start_pos_gap = read_start_pos - ref_start_pos
			end_pos_gap = read_end_pos - ref_end_pos
#			print("start_pos_gap: ", start_pos_gap, sep="\t")
#			print("end_pos_gap: ", end_pos_gap, sep="\t")
			gap_l.append(str(start_pos_gap))
			gap_l.append(str(end_pos_gap))
#		print("gap_l: ", gap_l, sep="\t")
		read_gap_dict[line_l[0]] = ",".join(gap_l)

	f4 = open(f_2)
	for line in f4:
		line = line.replace("\n","")
		line_l = line.split("\t")
#		print(line)
#		print(line_l[10], line_l[11], sep="\t")
##		if not line_l[11] == "*":
##			if int(line_l[10]) < int(line_l[11]):
##				continue
		print(line, read_gap_dict[line_l[-1]], sep="\t")

def merge_transcript3(f):

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

		
