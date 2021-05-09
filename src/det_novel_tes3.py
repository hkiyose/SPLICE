import sys
import re

f = open(sys.argv[1])
read_info_dict = {}
for line in f:
	line = line.replace("\n","")
	line_l = line.split("\t")
#	print(line)
	read_info_dict[line_l[0]] = ";".join(line_l)

f2 = open(sys.argv[2])
for line in f2:
	line = line.replace("\n","")
	line_l = line.split("\t")
#	print(line)
	read_l = line_l[-1].split(",")
	for i in range(len(read_l)):
		read_info_l = read_info_dict[read_l[i]].split(";")
#		print(read_l[i], read_info_l, sep="\t")
		
		if line_l[3] == "*":
			print("1", "0.00", line_l[2], read_info_l[14], line_l[5], line_l[6], read_info_l[12], read_info_l[13], line_l[12], read_l[i], sep="\t")
		else:
			print("1", "0.00", line_l[2], line_l[3], line_l[5], line_l[6], read_info_l[12], read_info_l[13], line_l[12], read_l[i], sep="\t")
#	print()
	
