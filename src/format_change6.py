import sys
import re

sample_count = 0
eva_known_full_dict = {}
transcript_read_count_dict = {}
f = open(sys.argv[1])
for line in f:
	line = line.replace("\n","")
	line_l = line.split("\t")
#	print(line)
	if line.startswith("CodingType"):
#		print(line)
#		print(len(line_l))
#		print(int((len(line_l) - 7)/2))
		sample_count = int((len(line_l) - 7)/2)
		print("#gene", "transcript", "known/novel", "coding/non-coding", "transcript_length", "novel_info", "\t".join(line_l[7:7+sample_count]), sep="\t")
		continue
	if line_l[2] == "SINGLE" or line_l[3] == "chrM":
		continue
#	print(line)
	if line_l[1] == "KNOWN":
		if line_l[2] == "FULL":
			eva_known_full_dict[line_l[5]] = 0
#			print("eva_known_full_dict: ", eva_known_full_dict, sep="\t")
	
#		print(line_l[0], line_l[1], "-", line_l[4], line_l[5], "-", line_l[7:7+sample_count], sep="\t")
		transcript_info = line_l[0] +";"+ line_l[1] + ";-;"+ line_l[4] +";"+ line_l[5] + ";-"
	else:
#		print(line_l[0], line_l[1], line_l[2], line_l[4], "-", line_l[6], line_l[7:7+sample_count], sep="\t")
		transcript_info = line_l[0] +";"+ line_l[1] +";"+ line_l[2] +";"+ line_l[4]  +";-;"+ line_l[6]
#	print("transcript_info: ", transcript_info, sep="\t")
	
	if not transcript_info in transcript_read_count_dict:
		transcript_read_count_dict[transcript_info] = line_l[7:7+sample_count]
#		print("transcript_read_count_dict: ", transcript_read_count_dict, sep="\t")
	else:
#		print("check")
#		print(transcript_read_count_dict[transcript_info])
		for i in range(0,sample_count):
#			print(i, line_l[7+i], transcript_read_count_dict[transcript_info][i], sep="\t")
			transcript_read_count_dict[transcript_info][i] = int(transcript_read_count_dict[transcript_info][i]) + int(line_l[7+i])
#		print(transcript_read_count_dict[transcript_info])			
#	print()	
	
for k,v in sorted(transcript_read_count_dict.items(), key=lambda x:int(x[1][0]), reverse=True):
	k_l = k.split(";")	
#	print(k,v,sep="\t")
#	print(k_l)
	v_s = []
	for i in range(len(v)):
		v_s.append(str(v[i]))
	eva_coding = "non-coding"
	if k_l[0] == "CODING":
		eva_coding = "coding"
#	print("eva_coding: ", eva_coding, sep="\t")
	if k_l[1] == "KNOWN":
#		print("KNOWN")
		eva_full = "partial"
		if k_l[4] in eva_known_full_dict:
			eva_full = "full"
#		print("eva_full: ", eva_full, sep="\t")
		print(k_l[3], k_l[4], "known", eva_coding, eva_full, "-", "\t".join(v_s), sep="\t")
	else:
#		print("NOVEL")
		eva_novel = ""
		if k_l[2] == "UNANNOTATED":
			eva_novel = "novel_exon"
		elif k_l[2] == "LENGTH":
			eva_novel = "novel_exon_length"
		elif k_l[2] == "COMBINATION":
			eva_novel = "novel_exon_combination"
#		print("eva_novel: ", eva_novel, sep="\t")
		print(k_l[3], k_l[4], eva_novel, eva_coding, "-", k_l[5], "\t".join(v_s), sep="\t")
#	print()
