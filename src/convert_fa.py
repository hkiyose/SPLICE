import sys
import re

seq = ""
f = open(sys.argv[1])
for line in f:
	line = line.replace("\n","")
	line_l = line.split("\t")
#	print(line)
	if line.startswith(">"):
		if not seq == "":
			print(seq)
			seq = ""
		print(line)
	else:
		seq += line
#		print(seq)
