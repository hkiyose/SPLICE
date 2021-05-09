import sys
import re
	
transcript_read_count_dict = {}
f = open(sys.argv[1])
for line in f:
	line = line.replace("\n","")
	line_l = line.split("\t")
#	print(line)
	k = "/".join(line_l[0:6])
#	print("k: ", k, sep="\t")

	if not k in transcript_read_count_dict:
		transcript_read_count_dict[k] = []
		for reads in line_l[6:]:
			transcript_read_count_dict[k].append(int(float(reads)))
#		print("transcript_read_count_dict[k]: ", transcript_read_count_dict[k], sep="\t")
	else:
		for i in range(len(line_l[6:])):
			transcript_read_count_dict[k][i] += int(float(line_l[i+6]))
#		print("MERGE")
#		print("transcript_read_count_dict[k]: ", transcript_read_count_dict[k], sep="\t")
#	print()

for k,v in transcript_read_count_dict.items():
	v_s = []
	for i in range(len(v)):
		v_s.append(str(v[i]))
	k_l = k.split("/")
	print("\t".join(k_l[0:5]), "/".join(k_l[5:]), "\t".join(v_s), sep="\t")	
