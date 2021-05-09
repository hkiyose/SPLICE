import sys
import re

splicing_variant_name_dict = {}
liver_count_dict = {}
cancer_count_dict = {}
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
	
	count = 0
	liver_count_dict[convert_sample_name] = []
	cancer_count_dict[convert_sample_name] = []
	for sample in line_l[10:]:
		count += 1
		if count % 2 == 1:
			liver_count_dict[convert_sample_name].append(sample)
		else:
			cancer_count_dict[convert_sample_name].append(sample)
	
	if line.startswith("CodingType"):
		continue
	else:
		print(convert_sample_name, line, sep="\t")
	
#	if line.startswith("CodingType"):	
#		print("TranscriptID", "\t".join(cancer_count_dict[convert_sample_name]), "\t".join(liver_count_dict[convert_sample_name]) ,sep="\t")
#	else:
#		print(convert_sample_name, "\t".join(cancer_count_dict[convert_sample_name]), "\t".join(liver_count_dict[convert_sample_name]) ,sep="\t")
#	print()
