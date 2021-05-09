import sys
import re

transcript_support_dict = {}
f = open(sys.argv[1])
for line in f:
	line = line.replace("\n","")
	line_l = line.split("\t")
#	print(line)
	if line.startswith("Variant_type"):
		print(line)
		continue
	info_l = line_l[0].split(";")
	if info_l[1] == "novel":
		continue
	gene_l = info_l[1].split("/")
	if len(gene_l) > 2:
		continue
		
	total_reads = 0 
	for i in range(len(line_l)):
		if i == 0:
			continue
		total_reads += int(line_l[i])
#	print(line)
#	print(total_reads)
	transcript_support_dict[line] = total_reads
#	print()
	
for k,v in sorted(transcript_support_dict.items(), key=lambda x:x[1], reverse=True):
	print(k)


