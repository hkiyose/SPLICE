import sys
import re

def del_term_index(l):
	if l[-1] == "":
		del l[-1]
	return l

gene_exon_num_pos_dict = {}
f = open(sys.argv[1])
for line in f:
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
f2 = open(sys.argv[2])
for line in f2:
	line = line.replace("\n","")
	line_l = line.split("\t")
	read_dict[line_l[-1]] = 0

read_gap_dict = {}
f3 = open(sys.argv[3])
for line in f3:
	line = line.replace("\n","")
	line_l = line.split("\t")
	
	if not line_l[0] in read_dict:
		continue	
	
#	print(line)
#	print(line_l[11], line_l[12], line_l[13], sep="\t")
#	print(line_l[6], line_l[7], sep="\t")

	exon_num_l = line_l[12].split(",")
	eva_exon_len_l = line_l[13].split(",")
	start_pos_l = line_l[6].split(",")
	end_pos_l = line_l[7].split(",")
#	print("exon_num_l: ", exon_num_l, "eva_exon_len_l: ", eva_exon_len_l, sep="\t")
#	print("start_pos_l: ", start_pos_l, "end_pos_l: ", end_pos_l, sep="\t")
	
	gap_l = []
	for i in range(len(exon_num_l)):
#		if i == 0 or i ==  len(exon_num_l)-1:
#			continue
		k = line_l[11] +"/"+ exon_num_l[i]
		
		ref_start_pos = int(gene_exon_num_pos_dict[k][0]) + 1
		ref_end_pos = int(gene_exon_num_pos_dict[k][1])
		read_start_pos = int(start_pos_l[i])
		read_end_pos = int(end_pos_l[i])
#		print("ref_start_pos: ", ref_start_pos, "ref_end_pos: ", ref_end_pos, "read_start_pos: ", read_start_pos, "read_end_pos: ", read_end_pos, sep="\t")
		
		start_pos_gap = read_start_pos - ref_start_pos
		end_pos_gap = read_end_pos - ref_end_pos
#		print("start_pos_gap: ", start_pos_gap, sep="\t")
#		print("end_pos_gap: ", end_pos_gap, sep="\t")
		gap_l.append(str(start_pos_gap))
		gap_l.append(str(end_pos_gap))
#	print("gap_l: ", gap_l, sep="\t")
	read_gap_dict[line_l[0]] = ",".join(gap_l)
	

f4 = open(sys.argv[2])
for line in f4:
	line = line.replace("\n","")
	line_l = line.split("\t")
#	print(line)
#	print(line_l[10], line_l[11], sep="\t")
##	if not line_l[11] == "*":
##		if int(line_l[10]) < int(line_l[11]):
##			continue
	print(line, read_gap_dict[line_l[-1]], sep="\t")

