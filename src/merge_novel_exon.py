import sys
import os
import re
from glob import glob
import math

dir_name = str(sys.argv[1])
file_l = []
for file in glob(dir_name + '/*.novel_exon'):
#	print(file)
	file_l.append(file)

variant_read_count_dict = {}#variant_read_count_dict[k] = [cancer_read_count, liver_read_count]
variant_error_rate_dict = {}
read_count_l = [0]*len(file_l)
sample_name_l = []
file_count = -1
for file in sorted(file_l):
	file_count += 1
	f = open(file)
	basename = os.path.basename(file)
#	print("basename: ", basename, sep="\t")
	basename_l = basename.split("_")
#	print(basename_l[0], basename_l[1], sep="\t")	
	
	sample_name = basename_l[0]
#	print("sample_name: ", sample_name, sep="\t")
	sample_name_l.append(sample_name)
	
	for line in f:
		line = line.replace("\n", "")
		line_l = line.split("\t")	
		read_count_l[file_count] += int(line_l[0])
		if line_l[0] == "1":
			continue
#		print(line)	

		k = line_l[1] +";"+ line_l[2] +";"+ line_l[3] +";"+ line_l[4]
#		print(k)
		
		if not k in variant_read_count_dict:
			variant_read_count_dict[k] = ["0"]*len(file_l)
			
#			print(variant_read_count_dict[k])
		
		variant_read_count_dict[k][file_count] = line_l[0]
#		print("sample_name_l: ", sample_name_l, sep="\t")
#		print("variant_read_count_dict[k]: ", variant_read_count_dict[k], sep="\t")		
#		print()

print("Variant_type","\t".join(sample_name_l),sep="\t")
for k,v in variant_read_count_dict.items():
	print(k,"\t".join(v),sep="\t")


