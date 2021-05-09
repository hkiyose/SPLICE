import sys
import os
import re
from glob import glob
import math

dir_name = str(sys.argv[1])
file_l = []
for file in glob(dir_name + '/*.annot'):
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
#		print(line)	
		read_count_l[file_count] += int(line_l[0])
		k = line_l[2] +"/"+ line_l[3] +"/"+ line_l[4] +"/"+ line_l[5] +"/"+ line_l[6] +"/"+ line_l[7] +"/"+ line_l[10]	
		if not k in variant_read_count_dict:
			variant_read_count_dict[k] = ["0"]*len(file_l)
			variant_error_rate_dict[k] = ["0"]*len(file_l)	
#			print(variant_read_count_dict[k])
#			print(variant_read_rate_dict[k])
		
		variant_read_count_dict[k][file_count] = line_l[0]
		variant_error_rate_dict[k][file_count] = line_l[8]
#		print("sample_name_l: ", sample_name_l, sep="\t")
#		print("variant_read_count_dict: ", variant_read_count_dict, sep="\t")		
#		print("variant_read_rate_dict: ", variant_read_rate_dict, sep="\t")
#		print()
#	print()
			
sample_name_read_rate_l = []
sample_name_error_rate_l = []
sample_name_fc_l = []
sample_name_log2_fc_l = []
for i in range(len(sample_name_l)):
	sample_name_read_rate_l.append(sample_name_l[i] + "_read/1Mread")
	sample_name_error_rate_l.append(sample_name_l[i] + "_error_rate")
	
	sample_name_l_l = sample_name_l[i].split("_")
	if not sample_name_l_l[0] in sample_name_fc_l:
		if not sample_name_l_l[0] + "_FC(Cancer/Liver)" in sample_name_fc_l:
			sample_name_fc_l.append(sample_name_l_l[0] + "_FC(Cancer/Liver)")
			sample_name_log2_fc_l.append(sample_name_l_l[0] + "_log2(FC(Cancer/Liver))")

print("Log2_FC_ave" , "variant_info", "\t".join(sample_name_l), "\t".join(sample_name_read_rate_l), "\t".join(sample_name_fc_l), "\t".join(sample_name_log2_fc_l), "\t".join(sample_name_error_rate_l), sep="\t")

#print(read_count_l)

for k,v in variant_read_count_dict.items():
#	print(len(read_count_l), read_count_l, sep="\t")
#	print(len(v), v, sep="\t")
	read_rate_l = [""]*len(file_l)
	for i in range(len(read_count_l)):
		if not v[i] == "0":
			read_rate_l[i] = str(int(v[i])/int(read_count_l[i])*1000000)
		else:
			read_rate_l[i] = "0"
#	print(len(read_rate_l), read_rate_l, sep="\t")

#	print(str(log2_fc_ave), k, "\t".join(v), "\t".join(read_rate_l), "\t".join(fc_l), "\t".join(log2_fc_l), "\t".join(variant_error_rate_dict[k]), sep="\t")
	print("-", k, "\t".join(v), "\t".join(read_rate_l), "-", "-", "\t".join(variant_error_rate_dict[k]), sep="\t")
#	print()
