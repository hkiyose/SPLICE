import sys
import re

f = open(sys.argv[1])
id_match_num_dict = {}
for line in f:
	line = line.replace("\n","")
	line_l = line.split("\t")
#	print(line)
	if not line_l[1] in id_match_num_dict:
		id_match_num_dict[line_l[1]] = line_l[0]
	else:
		if int(line_l[0]) > int(id_match_num_dict[line_l[1]]):#Select the one with the highest number of matches.
#			print("high")
			id_match_num_dict[line_l[1]] = line_l[0]
#	print(id_match_num_dict[line_l[1]])
#	print()
			
#for k,v in id_match_num_dict.items():
#	print(k,v,sep="\t")

f2 = open(sys.argv[2])
for line in f2:
	line = line.replace("\n","")
	line_l = line.split("\t")
#	print(line)	
	if line_l[0] in id_match_num_dict:
		print(line, id_match_num_dict[line_l[0]], sep="\t")
	else:
		print(line, "*", sep="\t")
	
