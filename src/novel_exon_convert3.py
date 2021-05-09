import sys
import re

def del_emp(l):
	l2 = []
	for i in range(len(l)):
		if not l[i] == "":
			l2.append(l[i])
	return l2

gene_pos_dict = {}
f = open(sys.argv[1])
for line in f:
	line = line.replace("\n","")
	line_l = line.split("\t")
	start_pos_l = del_emp(line_l[9].split(","))
	end_pos_l = del_emp(line_l[10].split(","))
	exon_num_l = line_l[-1].split(",")
	for i in range(len(start_pos_l)):
		k = line_l[12] +"/"+ exon_num_l[i]
		if not k in gene_pos_dict:
			gene_pos_dict[k] = [start_pos_l[i], end_pos_l[i]]
	
novel_transcript_eva_dict = {}
novel_transcript_convert_dict = {}
f2 = open(sys.argv[2])
for line in f2:
	line = line.replace("\n","")
	line_l = line.split("\t")
#	print(line)
	if line.startswith("TranscriptID"):
		continue
	gene_info = line_l[0].split("/")
	novel_pos_l = line_l[2].split("/")
#	print("novel_pos_l: ", novel_pos_l, sep="\t")
	novel_pos_l2 = novel_pos_l[1].split(",")
	exon_num_l = novel_pos_l[0].split(",")
	if "alt" in line_l[2] or "chrUn" in line_l[2] or "random" in line_l[2]:
#		print("chrUn")
		novel_transcript_eva_dict[gene_info[2] +"/"+ line_l[1]] = False
		continue
#	print("novel_pos_l2: ", novel_pos_l2, sep="\t")	
	start_l = []
	end_l = []
	novel_len_l = []
	for i in range(len(novel_pos_l2)):
#		print(novel_pos_l2[i])	
		if "chr" in novel_pos_l2[i]:
			pos_l = novel_pos_l2[i].split("_")
#			print("pos_l: ", pos_l, sep="\t")
			start_l.append(pos_l[1])
			end_l.append(pos_l[2])
			novel_len_l.append(int(pos_l[2]) - int(pos_l[1]) + 1)
		else:
			k = gene_info[2] + "/" + exon_num_l[i]
#			print("gene_pos_dict[k]: ", gene_pos_dict[k], sep="\t")
			start_l.append(gene_pos_dict[k][0])
			end_l.append(gene_pos_dict[k][1])
#	print("start_l: ", start_l, sep="\t")
#	print("end_l: ", end_l, sep="\t")
#	print("novel_len_l: ", novel_len_l, sep="\t")
	
	eva_novel_len = True
	for i in range(len(novel_len_l)):
		if int(novel_len_l[i]) < int(sys.argv[3]):
			eva_novel_len = False
#	print("eva_novel_len: ", eva_novel_len, sep="\t")
	
	eva_intron_len = True
	for i in range(len(start_l)):
		if i == 0:
			continue
		intron_len = int(end_l[i]) - int(start_l[i-1]) + 1
#		print("intron_len: ", intron_len, sep="\t")
		if intron_len < int(sys.argv[4]):
			eva_intron_len = False	
#	print("eva_intron_len: ", eva_intron_len, sep="\t")
		
	if eva_novel_len == False or eva_intron_len == False:
#		print("FilterOut")
		novel_transcript_eva_dict[gene_info[2] +"/"+ line_l[1]] = False
#		print()
		continue
	novel_transcript_eva_dict[gene_info[2] +"/"+ line_l[1]] = True
	novel_transcript_convert_dict[gene_info[2] +"/"+ line_l[1]] = line_l[2]
#	print()
	
#for k,v in novel_transcript_eva_dict.items():
#	print(k,v,sep="\t")#AFM/2,4,5,novel/s,k,k,chr4_73484454_73485677/*,*/*      True
	
f5 = open(sys.argv[5])
for line in f5:
	line = line.replace("\n","")
	line_l = line.split("\t")
#	print(line)
	if line_l[2] == "UNANNOTATED":
#		print(line)
#		print(line_l[3] +"/"+ line_l[5])
#		print(novel_transcript_eva_dict[line_l[3] +"/"+ line_l[5]])
		if novel_transcript_eva_dict[line_l[3] +"/"+ line_l[5]] == True:
#			print(novel_transcript_convert_dict[line_l[3] +"/"+ line_l[5]])
			print("\t".join(line_l[0:5]), novel_transcript_convert_dict[line_l[3] +"/"+ line_l[5]], "\t".join(line_l[6:]), sep="\t")
#		print()
#	print()
	else:
		print(line)



