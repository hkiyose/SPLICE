import sys
import re

def del_terminal_index(l):
	if l[-1] == "":
		del l[-1]
	return l

f = open(sys.argv[1])
chr_seq_dict = {}
for line in f:
	line = line.replace("\n", "")
	line_l = line.split()
	if line.startswith(">"):
		if line.startswith(">NC_"):
			chr2 = line_l[0][1:]
		else:
			chr2 = line[1:]
		chr_seq_dict[chr2] = ""
	else:
		chr_seq_dict[chr2] = line

		f2 = open(sys.argv[2])
		seq = ""
		for line in f2:
			line = line.replace("\n", "")
			line_l = line.split("\t")
			start_pos_l = del_terminal_index(line_l[9].split(","))
			end_pos_l = del_terminal_index(line_l[10].split(","))
			if line_l[2] in chr_seq_dict:
				for i in range(len(start_pos_l)):
					seq += chr_seq_dict[line_l[2]][int(start_pos_l[i]):int(end_pos_l[i])]
				if not seq == "":
					print(">" + line_l[1] + "/" + line_l[12])
					print(seq)
			seq = ""
		chr_seq_dict = {}
