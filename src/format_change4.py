import sys
import re

f = open(sys.argv[1])
for line in f:
	line = line.replace("\n","")
	line_l = line.split("\t")
#	print(line)
	if line.startswith("CodingType"):
		print(line)
		continue
#	print(line_l[8], line_l[9])
	if float(line_l[8]) < float(sys.argv[2]) or int(line_l[9]) < int(sys.argv[3]):
		continue
	print(line)
#	print()
