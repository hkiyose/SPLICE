import sys
import re

def del_term_index(l):
	if l[-1] == "":
		del l[-1]
	return l

gene_pos_dict = {}#gene_pos_dict[gene/exon_num] = [start_pos, end_pos]
f = open(sys.argv[1])
for line in f:
	line = line.replace("\n", "")
	line_l = line.split("\t")
#	print(line)
	start_pos_l = del_term_index(line_l[9].split(","))
	end_pos_l = del_term_index(line_l[10].split(","))
	exon_num_l = line_l[-1].split(",")
#	print("start_pos_l: ", start_pos_l, sep="\t")
#	print("end_pos_l: ", end_pos_l, sep="\t")
#	print("exon_num_l: ", exon_num_l, sep="\t")
	for i in range(len(start_pos_l)):
		k = line_l[12] +"/"+ exon_num_l[i]
		if not k in gene_pos_dict:
			gene_pos_dict[k] = [start_pos_l[i], end_pos_l[i]]

#for k,v in gene_pos_dict.items():
#	print(k,v,sep="\t")	

pre_gene = ""
f2 = open(sys.argv[2])
for line in f2:
	line = line.replace("\n", "")
	line_l = line.split("\t")
#	print(line)
	if "/" in line_l[3]:	
#		print()
		continue
	if "Log2_FC_ave" in line:
		continue
#	print(line)
	info_l = line_l[5].split("/")
	exon_num_l = info_l[0].split(",")
#	print("info_l: ", info_l, sep="\t")
	eva_exon_len_l = info_l[1].split(",")
#	print("info_l: ", info_l, sep="\t")
#	print("exon_num_l: ", exon_num_l, sep="\t")
#	print("eva_exon_len_l: ", eva_exon_len_l, sep="\t")

	start_pos_l2 = []
	end_pos_l2 = []
	
	transcript_length = 0
	for i in range(len(exon_num_l)):
		if exon_num_l[i] == "novel":
			continue
#		print(exon_num_l[i], gene_pos_dict[line_l[3]+"/"+exon_num_l[i]], sep="\t")
		transcript_length += int(gene_pos_dict[line_l[3]+"/"+exon_num_l[i]][1]) - int(gene_pos_dict[line_l[3]+"/"+exon_num_l[i]][0]) + 1
		start_pos_l2.append(int(gene_pos_dict[line_l[3]+"/"+exon_num_l[i]][0]))
		end_pos_l2.append(int(gene_pos_dict[line_l[3]+"/"+exon_num_l[i]][1]))
	
	for i in range(len(eva_exon_len_l)):
		if not "_" in eva_exon_len_l[i]:
			continue
#		print(eva_exon_len_l[i])
		eva_exon_len_l_l = eva_exon_len_l[i].split("_")
#		print(eva_exon_len_l_l[-2], eva_exon_len_l_l[-1], sep="\t")
		transcript_length += int(eva_exon_len_l_l[-1]) - int(eva_exon_len_l_l[-2]) + 1
		start_pos_l2.append(int(eva_exon_len_l_l[-2]))
		end_pos_l2.append(int(eva_exon_len_l_l[-1]))	
	
#	print("start_pos_l2: ", sorted(start_pos_l2), sep="\t")	
#	print("end_pos_l2: ", sorted(end_pos_l2), sep="\t")
#	print("transcript_length: ", transcript_length, sep="\t")
	
	start_pos_l3 = []
	end_pos_l3 = []
	for pos in sorted(start_pos_l2):
		start_pos_l3.append(str(pos))	
	for pos in sorted(end_pos_l2):
		end_pos_l3.append(str(pos))
	
	print("\t".join(line_l[0:6]), transcript_length, ",".join(start_pos_l3), ",".join(end_pos_l3), "\t".join(line_l[6:]), sep="\t")
		
#	print()
