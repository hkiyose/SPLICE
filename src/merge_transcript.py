import sys
import re

def del_terminal_index(l):
	if l[-1] == "":
		del l[-1]
	return l

f = open(sys.argv[1])
readID_merge_pos_dict = {}
readID_junction_match_num_dict = {}
readID_genome_transcriptome_match_num_dict = {}
for line in f:
	line = line.replace("\n", "")
	line_l = line.split("\t")
	start_pos_l = line_l[6].split(",")
	end_pos_l = line_l[7].split(",")
#	print(line)
	readID_merge_pos_dict[line_l[0]] = [end_pos_l[0], start_pos_l[-1]]
	readID_junction_match_num_dict[line_l[0]] = line_l[16]
	readID_genome_transcriptome_match_num_dict[line_l[0]] = [line_l[17],line_l[18]]
#       print(readID_merge_pos_dict[line_l[0]])
	
f2 = open(sys.argv[2])
for line in f2:
	line = line.replace("\n", "")
	line_l = line.split("\t")
	exon_num_l = line_l[5].split(",")	
	eva_exonlen_l = del_terminal_index(line_l[6].split(","))
	readID_l = line_l[9].split(",")
##	print(line)
	
	for i in range(len(readID_l)):
#		print("1", "0.0", line_l[2], line_l[3], line_l[4], gene_exonNum_pos_dict[line_l[2] +"," + exon_num_l[0]][1], gene_exonNum_pos_dict[line_l[2] +","+ exon_num_l[-1]][0], "*,"+",".join(exon_num_l[1:-1])+",*", "*,"+",".join(eva_exonlen_l[1:-1])+",*", line_l[9], sep="\t" )
		print("1", "0.0", line_l[2], line_l[3], line_l[4], readID_merge_pos_dict[readID_l[i]][0], readID_merge_pos_dict[readID_l[i]][1], "*,"+",".join(exon_num_l[1:-1])+",*", "*,"+",".join(eva_exonlen_l[1:-1])+",*", readID_junction_match_num_dict[readID_l[i]], readID_genome_transcriptome_match_num_dict[readID_l[i]][0], readID_genome_transcriptome_match_num_dict[readID_l[i]][1], readID_l[i], sep="\t" )
		
##	print()
