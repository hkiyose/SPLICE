import sys
import re

f = open(sys.argv[1])
for line in f:
	line = line.replace("\n","")
	line_l = line.split("\t")
#	print(line)

	chr_l = line_l[3].split("/")
	breakpoint_l = line_l[4].split("/")
	breakpoint_gap = abs(int(breakpoint_l[0]) - int(breakpoint_l[1]))	
#	print(line)
#	print("breakpoint_gap: ", breakpoint_gap, sep="\t")

	if chr_l[0] == chr_l[1]:
		if breakpoint_gap <= int(sys.argv[2]):
#			print(">>>close")
#			print()
			continue
		
	print(line)
	
#	print()
	
		







