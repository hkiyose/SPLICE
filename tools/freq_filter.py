import sys
import re

def freq_filter(f, int2, int3):

	gene_read_count_dict = {}
	sample_count = 0
	f1 = open(f)
	for line in f1:
		line = line.replace("\n","")
		line_l = line.split("\t")
	#	print(line)
		if line.startswith("Log2_FC_ave"):
			print(line)
			for i in range(len(line_l)):
				if "_read/1Mread" in line_l[i]:
					sample_count = i
					break
			continue
	#	print("sample_count: ", sample_count, sep="\t")
	#	print(line)
		read_l = line_l[2:sample_count]
	#	print("read_l: ", read_l, sep="\t")
		total_read_count = 0
		for i in range(len(read_l)):
			total_read_count += int(float(read_l[i]))
	#	print("total_read_count: ", total_read_count, sep="\t")
		info_l = line_l[1].split("/")
		if info_l[0] in gene_read_count_dict:
			gene_read_count_dict[info_l[0]] += total_read_count
		else:
			gene_read_count_dict[info_l[0]] = total_read_count

	#for k,v in gene_read_count_dict.items():
	#	print(k,v,sep="\t")

	f_2 = open(f)
	for line in f_2:
		line = line.replace("\n","")
		line_l = line.split("\t")
	#	print(line)
		if line.startswith("Log2_FC_ave"):
			continue
		info_l = line_l[1].split("/")
		read_l = line_l[2:sample_count]
	#	print("read_l: ", read_l, sep="\t")
		total_read_count = 0
		max_read_count = 0	
		for i in range(len(read_l)):
			total_read_count += int(float(read_l[i]))
			if int(float(read_l[i])) > max_read_count:
				max_read_count = int(float(read_l[i]))
	#	print("total_read_count: ", total_read_count, "splicing_variant_freq: ", round(total_read_count/gene_read_count_dict[info_l[0]]*100, 3), sep="\t")
	#	print("max_read_count: ", max_read_count, sep="\t")
		if info_l[1] == "NOVEL":
			if max_read_count >= int(int2) and total_read_count/gene_read_count_dict[info_l[0]]*100 >= int(int3):
				print(line)
		else:
			print(line)

			

	
def freq_filter2(f, f_2):

	gene_read_count_dict = {}
	sample_count = 0
	f1 = open(f)
	for line in f1:
		line = line.replace("\n","")
		line_l = line.split("\t")
		if line.startswith("Log2_FC_ave"):
			print(line)
			for i in range(len(line_l)):
				if "_read/1Mread" in line_l[i]:
					sample_count = i
			continue
		read_l = line_l[2:sample_count]
		total_read_count = 0
		for i in range(len(read_l)):
			total_read_count += int(float(read_l[i]))
		info_l = line_l[1].split("/")	
		if info_l[0] in gene_read_count_dict:
			gene_read_count_dict[info_l[0]] += total_read_count
		else:
			gene_read_count_dict[info_l[0]] = total_read_count

	#for k,v in gene_read_count_dict.items():
	#	print(k,v,sep="\t")

	f2 = open(f_2)
	for line in f2:
		line = line.replace("\n","")
		line_l = line.split("\t")
	#	print(line)
		if line_l[0] == "NOVEL" or "/" in line_l[0]:
			print(line)
			continue
		read_l = line_l[3:]
	#	print(read_l)
		total_read_count = 0
		max_read_count = 0	
		for i in range(len(read_l)):
			total_read_count += int(float(read_l[i]))
			if int(float(read_l[i])) > max_read_count:
				max_read_count = int(float(read_l[i]))
	#	print("total_read_count: ", total_read_count, round(total_read_count/gene_read_count_dict[line_l[3]]*100, 3), sep="\t")
	#	print("max_read_count: ", max_read_count, sep="\t")
		if line_l[0] in gene_read_count_dict:
			if max_read_count >= 3 and total_read_count/gene_read_count_dict[line_l[0]]*100 >= 1:
				print(line)
		
