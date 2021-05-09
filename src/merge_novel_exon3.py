import sys
import re

merged_pos = []
f = open(sys.argv[1])
num_lines = sum(1 for line in open(sys.argv[1]))
#print(num_lines)
line_count = 0
for line in f:
	line_count += 1
	line = line.replace("\n","")
	line_l = line.split("\t")
#	print("line: ", line, sep="\t")
	
	if len(merged_pos) == 0:
		merged_pos = [line_l[0], line_l[1], line_l[2]]
#		print("merged_pos: " ,merged_pos ,sep="\t")
#		print()
		continue
	
	if int(line_l[1]) > int(merged_pos[2]) or line_l[0] != merged_pos[0]:
#		print("!")
		print("\t".join(merged_pos))
		merged_pos = [line_l[0], line_l[1], line_l[2]]
#		print("merged_pos: " , merged_pos, sep="\t")
		if num_lines == line_count:
			print("\t".join(merged_pos))
	else:
		novel_pos_in_merged_set = set([x for x in range(int(merged_pos[1]), int(merged_pos[2])+1)])
		novel_pos_in_line_set = set([x for x in range(int(line_l[1]), int(line_l[2])+1)])
		and_set = novel_pos_in_merged_set & novel_pos_in_line_set
#		print("novel_pos_in_merged_set: ", len(novel_pos_in_merged_set), novel_pos_in_merged_set, sep="\t")
#		print("novel_pos_in_line_set: ", len(novel_pos_in_line_set), novel_pos_in_line_set, sep="\t")
#		print("and_set_count: ", len(and_set), and_set, sep="\t")
		if len(and_set) > 0:	
#			print("overlap")
			if int(line_l[2]) > int(merged_pos[2]):
				merged_pos[2] = line_l[2]
#				print("merged_pos: " , merged_pos, sep="\t")
		if num_lines == line_count:
			print("\t".join(merged_pos))
#	print()	
	

