import sys
import re

def del_term_index(l):
	if l[-1] == "":
		del l[-1]
	return l

def del_empty_index(l):
	l2 = []
	for i in range(len(l)):
		if not l[i] == "" and not l[i] == "*":
			l2.append(l[i])
	return l2

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

f2 = open(sys.argv[2])
for line in f2:
	line = line.replace("\n","")
	line_l = line.split("\t")
	if not line_l[1] == "NOVEL":
		continue
	if not line_l[2] == str(sys.argv[3]):
		continue
#	print(line)
	info_l = line_l[5].split("/")
	exon_num_l = del_empty_index(info_l[0].split(","))
	merged_exon_pos = info_l[2].split(",")
	gap_l = del_empty_index(info_l[3].split(","))
#	print("info_l: ", info_l, sep="\t")
#	print("exon_num_l: ", exon_num_l, sep="\t")
#	print("merged_exon_pos: ", merged_exon_pos, sep="\t")
#	print("gap_l: ", gap_l, sep="\t")
	
	pos_l = []	
	pos_l.append(merged_exon_pos[0])
	for i in range(len(exon_num_l)):
#		print(gene_pos_dict[line_l[3]+"/"+exon_num_l[i]])
		pos_l.append(gene_pos_dict[line_l[3]+"/"+exon_num_l[i]][0])
		pos_l.append(gene_pos_dict[line_l[3]+"/"+exon_num_l[i]][1])
	pos_l.append(merged_exon_pos[1])
	
	pos_l2 = []
	start_pos_l = []
	end_pos_l = []	
#	print("pos_l: ", pos_l, sep="\t")
	for i in range(len(pos_l)):
#		print("pos: ", pos_l[i], "gap: ", gap_l[i], sep="\t")
		if i == 0 or i == len(pos_l)-1:
#			print("nochange")
			pos_l2.append(pos_l[i])
		else:
			pos_l2.append(str(int(pos_l[i])+int(gap_l[i])))
			if i%2 == 1:
				start_pos_l.append(str(int(pos_l[i])+int(gap_l[i])))
			else:
				end_pos_l.append(str(int(pos_l[i])+int(gap_l[i])))
#	print("pos_l2: ", pos_l2, sep="\t")
#	print("start_pos_l: ", start_pos_l, sep="\t")
#	print("end_pos_l: ", end_pos_l, sep="\t")
	
	transcript_length = 0
	for i in range(len(start_pos_l)):
		transcript_length += int(end_pos_l[i]) - int(start_pos_l[i]) + 1
#	print("transcript_length: ", transcript_length, sep="\t")	
	print("\t".join(line_l[0:6]), transcript_length, ",".join(pos_l2), "\t".join(line_l[6:]), sep="\t")
#	print()
	
