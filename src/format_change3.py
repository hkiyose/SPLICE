import sys
import re

sample_count = 0
f = open(sys.argv[1])
for line in f:
	line = line.replace("\n","")
	line_l = line.split("\t")
	if line.startswith("Log2_FC_ave"):
		for i in range(len(line_l)):
			if "_read/1Mread" in line_l[i]:
				sample_count = i
		print("CodingType", "EvaKnown", "EvaLength", "Gene", "SplicingVariant", "VariantInfo", "TotalReadCount", "TotalGeneCount", "SplicingVariantFrequency(%)", "MaxReadCount", "\t".join(line_l[2:sample_count+1]), sep="\t")
		break
	break

sample_count = 0	
gene_read_count_dict = {}
f2 = open(sys.argv[2])
for line in f2:
	line = line.replace("\n","")
	line_l = line.split("\t")
#	print(line)
	if sample_count == 0:
		sample_count = int((len(line_l)-6)/2)
#		print("sample_count: ", sample_count)
	read_l = line_l[6:6+sample_count]
#	print(read_l)
	total_read_count = 0
	for i in range(len(read_l)):
#		print(read_l[i])
		total_read_count += int(read_l[i])
#	print("total_read_count: ", total_read_count)
	if line_l[3] in gene_read_count_dict:
		gene_read_count_dict[line_l[3]] += total_read_count
	else:
		gene_read_count_dict[line_l[3]] = total_read_count
#	print()
	
#for k,v in gene_read_count_dict.items():
#	print(k,v,sep="\t")

f3 = open(sys.argv[2])
for line in f3:
	line = line.replace("\n","")
	line_l = line.split("\t")
#	print(line)
	read_l = line_l[6:6+sample_count]
	total_read_count = 0
	max_read_count = 0	
	for i in range(len(read_l)):
		total_read_count += int(read_l[i])
		if int(read_l[i]) > max_read_count:
			max_read_count = int(read_l[i])
#	print("total_read_count: ", total_read_count, round(total_read_count/gene_read_count_dict[line_l[3]]*100, 3), sep="\t")
#	print("max_read_count: ", max_read_count, sep="\t")
	print("\t".join(line_l[0:6]), total_read_count, gene_read_count_dict[line_l[3]], round(total_read_count/gene_read_count_dict[line_l[3]]*100, 3), max_read_count, "\t".join(line_l[6:]), sep="\t")
#	print()
		
