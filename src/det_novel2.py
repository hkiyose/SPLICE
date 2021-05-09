import sys
import re
import align2

def del_emp(l):
	l2 = []
	for i in range(len(l)):
		if not l[i] == "":		
			l2.append(l[i])
	return l2
	
def get_primerSeq(query_seq):
	primer_seq = "AAGCAGTGGTATCAACGCAGAGTAC"#SMARTer II A Oligonucleotide
	rev_primer_seq = "GTACTCTGCGTTGATACCACTGCTT"
		
	res1 = align2.SW(primer_seq, query_seq, 1, -1, -1)
	res2 = align2.SW(rev_primer_seq, query_seq, 1, -1, -1)
	consensus1 = ""
	match1 = 0
	
	for i in range(len(res1[0])):
		if res1[0][i] == res1[1][i]:
			match1 += 1
			consensus1 += "|"
		else:
			consensus1 += "*"
	
#	print(res1[0])
#	print(consensus1)
#	print(res1[1])
	
	consensus2 = ""
	match2 = 0
	
	for i in range(len(res2[0])):
		if res2[0][i] == res2[1][i]:
			match2 += 1
			consensus2 += "|"
		else:
			consensus2 += "*"
	
#	print(res2[0])
#	print(consensus2)
#	print(res2[1])
	
#	print("match1: ", match1, match1/len(res1[0]), sep="\t")
#	print("match2: ", match2, match2/len(res2[0]), sep="\t")
		
	if match1 >= int(sys.argv[7]) and match1/len(res1[0]) >= float(sys.argv[8]) or match2 >= int(sys.argv[7]) and match2/len(res2[0]) >= float(sys.argv[8]):
#		print("MATCH!")
		return "match"
	else:
#		print("unMATCH!!")
		return "unmatch"


f2 = open(sys.argv[3])#cDNAへのmapping結果
high_map_rate_dict = {}
for line in f2:
	line = line.replace("\n","")
	line_l = line.split("\t")
#	print(line)
	if float(line_l[1]) >= float(sys.argv[4]):
		high_map_rate_dict[line_l[2]] = 0

f = open(sys.argv[1])
novel_dict = {}
novel_read_dict = {}
fusion_id_seq_dict = {}#190406
fusion_id_info_dict = {}#190406
for line in f:
	line = line.replace("\n","")
	line_l = line.split("\t")
	
#	print(line)
#	print(line_l[11])
	
#	if not line_l[0] in unmapped_id_dict:
#		continue
	
	if "novel" in line_l[11]:
#		print(line)
#		print(line_l[11])
		eva_mq = "high"
		mq_l = line_l[9].split("/")
		for i in range(len(mq_l)):
			if not int(mq_l[i]) >= int(sys.argv[2]):
				eva_mq = "low"
		if eva_mq == "low":#mqが低いもの
#			print("eva_mq: ", eva_mq)
			continue
		
		eva_high_map_rate_in_cdna = False
		eva_known_gene = False
		annot_gene_l = del_emp(line_l[11].split("/"))
		for i in range(len(annot_gene_l)):
			if annot_gene_l[i] != "novel":#Fusions with known genes
				eva_known_gene = True
		
		if eva_known_gene == False:
			if line_l[0] in high_map_rate_dict:
				eva_high_map_rate_in_cdna = True
		
		if eva_high_map_rate_in_cdna == True:#Only novel exon is annotated & mapped with high mapping rate in cDNA		
#			print("high_map_rate_in_cdna")
			continue
	
#		print(line)		

		chr_l = line_l[3].split("/")
		start_pos_l = re.split('[,/]', line_l[6])
		end_pos_l = re.split('[,/]', line_l[7])
		annot_gene_l = del_emp(line_l[11].split("/"))
		exon_num_l = del_emp(line_l[12].split("/"))
		annot_l = del_emp(re.split('[,/]', line_l[13]))
		
		chr_l2 = []
		start_pos_l2 = line_l[6].split("/")
		for i in range(len(start_pos_l2)):
			start_pos_l2_l = start_pos_l2[i].split(",")
			for j in range(len(start_pos_l2_l)):
				chr_l2.append(chr_l[i])
#		print("start_pos_l: ", start_pos_l, sep="\t")
#		print("end_pos_l: ", end_pos_l, sep="\t")
#		print("annot_l: ", annot_l, sep="\t")
#		print("annot_gene_l: ", annot_gene_l, sep="\t")
#		print("chr_l: ", chr_l, sep="\t")
#		print("chr_l2: ", chr_l2, sep="\t")
		
		annot_l2 = []
		if not len(annot_l) == len(start_pos_l):
			continue
		for i in range(len(start_pos_l)):
			if annot_l[i] == "novel":	
				k = chr_l2[i]+"_"+start_pos_l[i]+"_"+end_pos_l[i]
				annot_l2.append(k)

#				print(k)
##				if k in novel_exon_dict:
##					novel_exon_dict[k] += 1
##				else:
##					novel_exon_dict[k] = 1
			else:
				annot_l2.append(annot_l[i])
		
#		if len(annot_gene_l) >= 2: 
#			print(line)
#			print(line_l[0], "/".join(chr_l), "/".join(annot_gene_l), "/".join(exon_num_l), ",".join(annot_l2), sep="\t")
#			print("")
		
#			print("chr_l: ", chr_l, sep="\t")
#			print("annot_gene_l: ", annot_gene_l, sep="\t")
#			print("exon_num_l: ", exon_num_l, sep="\t")
#			print("annot_l2: ", annot_l2, sep="\t")
#			if "/" in line_l[4]:
#				print("FUSION?")
#			print()

		gene_exon_num_l = []
		for i in range(len(annot_gene_l)):
			gene_exon_num = annot_gene_l[i] +"/"+ exon_num_l[i]
			gene_exon_num_l.append(gene_exon_num)  
		
		gene_exon_num_l_sort = sorted(gene_exon_num_l)	
#		print("gene_exon_num_l_sort: ", gene_exon_num_l_sort, sep="\t")
		
		annot_gene_l2 = []
		exon_num_l2 = []
		for i in range(len(gene_exon_num_l_sort)):
			gene_exon_num_l_sort_l = gene_exon_num_l_sort[i].split("/")
			annot_gene_l2.append(gene_exon_num_l_sort_l[0])
			exon_num_l2.append(gene_exon_num_l_sort_l[1])

		k = "/".join(chr_l)+":"+"/".join(annot_gene_l)+":"+ "/".join(exon_num_l)+":"+",".join(annot_l2)
#		print("k: ", k, sep="\t")
		#chr19/chr19:novel/PLA2G4C:*/61:chr19_48108487_48108663,long
#		print("")
			
		eva_novel = False#If annot_l2 is empty or only known, it will be false and will not be added to novel_dict.
		for i in range(len(annot_l2)):
			if not annot_l2[i] == "known" and not annot_l2[i] == "long" and not annot_l2[i] == "short" and not annot_l2[i] == "slong":
				eva_novel = True
#			print("eva_novel: ", eva_novel, sep="\t")
		
#		print("eva_novel: ", eva_novel, sep="\t")
		
		if "/" in line_l[4]:#Sort the order of annot_gene_l, etc. (for merge)
#			print("fusion?")
#			print(line)
#			print(line_l[0], "/".join(chr_l), "/".join(annot_gene_l), "/".join(exon_num_l), ",".join(annot_l2), sep="\t")
			annot_gene_l_tmp = [annot_gene_l[0], annot_gene_l[-1]]
#			print([annot_gene_l[0], annot_gene_l[-1]], sorted(annot_gene_l_tmp), sep="\t")
			
			if [annot_gene_l[0], annot_gene_l[-1]] == sorted(annot_gene_l_tmp) and annot_gene_l[0] != annot_gene_l[-1]:
				k = "/".join(chr_l)+":"+"/".join(annot_gene_l)+":"+ "/".join(exon_num_l)+":"+",".join(annot_l2)
			elif [annot_gene_l[0], annot_gene_l[-1]] == sorted(annot_gene_l_tmp) and annot_gene_l[0] == annot_gene_l[-1]:
				chr_l_tmp = [chr_l[0], chr_l[-1]]
#				print(chr_l_tmp, sorted(chr_l_tmp), sep="\t")
				if chr_l_tmp == sorted(chr_l_tmp):	
					k = "/".join(chr_l)+":"+"/".join(annot_gene_l)+":"+ "/".join(exon_num_l)+":"+",".join(annot_l2)
				else:
					k = "/".join(reversed(chr_l))+":"+"/".join(reversed(annot_gene_l))+":"+"/".join(reversed(exon_num_l))+":"+",".join(reversed(annot_l2))
			else:
				k = "/".join(reversed(chr_l))+":"+"/".join(reversed(annot_gene_l))+":"+"/".join(reversed(exon_num_l))+":"+",".join(reversed(annot_l2))
#			print("k: ", k, sep="\t")
#			print()	
			
		else:
			k = "/".join(chr_l)+":"+"/".join(annot_gene_l)+":"+ "/".join(exon_num_l)+":"+",".join(annot_l2)
#			print("k: ", k, sep="\t")
#			print()
		
		if eva_novel == True:
			if k in novel_dict:
#				print(novel_dict[k])
				novel_dict[k] += 1
				novel_read_dict[k] += ","+line_l[0]
				if "/" in line_l[4]:#Check if there is a primer sequence at breakpoint in the next step.
#					print("fusion?")
					fusion_id_seq_dict[line_l[0]] = 0
					fusion_id_info_dict[line_l[0]] = line
			else:
				novel_dict[k] = 1
				novel_read_dict[k] = line_l[0]
#				print(novel_dict[k])
				if "/" in line_l[4]:
#					print("fusion?")
					fusion_id_seq_dict[line_l[0]] = 0
					fusion_id_info_dict[line_l[0]] = line
#			print()			


##for k,v in novel_exon_dict.items():
##	print(k,v,sep="\t")

#Get the candidate sequences that might be fusion, and check if there is a primer sequence at breakpoint.
line_num = 0
read_id = ""
f3 = open(sys.argv[6])
for line in f3:
	line = line.replace("\n","")
	line_l = line.split()
	line_num += 1
	if line_num % 4 == 1:
		if line_l[0][1:] in fusion_id_seq_dict:
#			print(line)
			read_id = line_l[0][1:]
		else:
			read_id = ""
	elif line_num % 4 == 2 and not read_id == "":
		fusion_id_seq_dict[read_id] = line
#		print(fusion_id_seq_dict[read_id]) 

artifact_fusion_dict = {}
for k,v in fusion_id_seq_dict.items():
#	print(k,v)
#	print(fusion_id_info_dict[k])
	info_l = fusion_id_info_dict[k].split("\t")

	start_pos_l = info_l[4].split("/")#read上
	start_pos_l2 = start_pos_l[1].split(",")
	breakpoint_pos_end = start_pos_l2[0]

	end_pos_l = info_l[5].split("/")
	end_pos_l2 = end_pos_l[0].split(",")
	breakpoint_pos_start = end_pos_l2[-1]
	
#	print("breakpoint_pos_start: ", breakpoint_pos_start, sep="\t")  
#	print("breakpoint_pos_end: ", breakpoint_pos_end, sep="\t")
	if int(breakpoint_pos_end) - int(breakpoint_pos_start) >= int(sys.argv[5]):
#		print("60bp以上")
#		print(v)
#		print(v[int(breakpoint_pos_start)+1:int(breakpoint_pos_end)])
		if get_primerSeq(v[int(breakpoint_pos_start)+1:int(breakpoint_pos_end)]) == "match":
#			print("artifact_fusion")
			artifact_fusion_dict[k] = 0
#	print()

for k,v in sorted(novel_dict.items(), key=lambda x: -x[1]):
#	print(k,v)
	k_l = k.split(":")
	read_l = novel_read_dict[k].split(",")
#	print(v, "\t".join(k_l), novel_read_dict[k], sep="\t")
	read_l_rm_artifact_fusion = ""
	for i in range(len(read_l)):
		if not read_l[i] in artifact_fusion_dict:
#			print("ARTIFACT", read_l[i], sep="\t")			
#		else:
			if read_l_rm_artifact_fusion == "":
				read_l_rm_artifact_fusion += read_l[i]	
			else:
				read_l_rm_artifact_fusion += "," + read_l[i]
	
	read_l_rm_artifact_fusion_l = read_l_rm_artifact_fusion.split(",")
#	print("read_l: ", len(read_l))
#	print("read_l_rm_artifact_fusion_l: ", len(read_l_rm_artifact_fusion_l))
	if not read_l_rm_artifact_fusion_l[0] == "":
		print(len(read_l_rm_artifact_fusion_l), "\t".join(k_l), read_l_rm_artifact_fusion, sep="\t")
#	else:	
#		print(v, "\t".join(k_l), novel_read_dict[k], sep="\t")
#		print("read_l: ", len(read_l))
#		print("read_l_rm_artifact_fusion_l: ", len(read_l_rm_artifact_fusion_l))
				
#			print(read_l[i])
#	print("")
