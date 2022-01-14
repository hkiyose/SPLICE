import sys
import re

splicing_variant_name_dict = {}
f = open(sys.argv[1])
for line in f:
	line = line.replace("\n","")
	line_l = line.split("\t")
#	print(line)
	if line.startswith("Gene"):
		continue
	sample_name = "/".join(line_l[1:5])
#	print(sample_name)
	if sample_name in splicing_variant_name_dict:
		splicing_variant_name_dict[sample_name] += 1
	else:
		splicing_variant_name_dict[sample_name] = 0
	convert_sample_name = sample_name +"/"+ str(splicing_variant_name_dict[sample_name])
#	print("convert_sample_name: " , convert_sample_name, sep="\t")
	if line.startswith("CodingType"):
		continue
	else:
		print(convert_sample_name, line, sep="\t")
	
