import sys
import re

chr_pos_d = {}
f = open(sys.argv[1])
for line in f:
	line = line.replace("\n","")
	line_l = line.split("\t")
#	print(line)
	v = [line_l[1], line_l[2]]
	if line_l[0] in chr_pos_d:
		chr_pos_d[line_l[0]].append(v)	
	else:
		chr_pos_d[line_l[0]] = []
		chr_pos_d[line_l[0]].append(v)

#for k,v in chr_pos_d.items():
#	print(k,v,sep="\t")

novel_pos_dict = {}
f2 = open(sys.argv[2])
for line in f2:
	line = line.replace("\n","")
	line_l = line.split("\t")
#	print(line)
	if line.startswith("Vari"):
		print(line)
		continue
	
	info_l = line_l[0].split(";")
#	print("info_l: ", info_l, sep="\t")
	
	novel_pos_l = info_l[-1].split(",")
#	print("novel_pos_l: ", novel_pos_l, sep="\t")

	converted_novel_pos_l = []
	
	for i in range(len(novel_pos_l)):
#		print(novel_pos_l[i])
		if not "_" in novel_pos_l[i]:
			converted_novel_pos_l.append(novel_pos_l[i])
		else:
#			print(novel_pos_l[i])
			pos_l = novel_pos_l[i].split("_")
#			print("pos_l: ", pos_l, sep="\t")	
			chr_s = ""
			if len(pos_l) == 5:	
				chr_s = "_".join(pos_l[0:3])
			elif len(pos_l) == 4:
				chr_s = "_".join(pos_l[0:2])
			else:
				chr_s = "_".join(pos_l[0:1])	
#			print("chr_s: ", chr_s,sep="\t")
	
			for j in range(len(chr_pos_d[chr_s])):
#				print(chr_pos_d[chr_s][j][0], chr_pos_d[chr_s][j][1], sep="\t")
				if int(chr_pos_d[chr_s][j][0]) <= int(pos_l[-2]) <= int(chr_pos_d[chr_s][j][1]):
#					print(chr_pos_d[chr_s][j][0], chr_pos_d[chr_s][j][1], sep="\t")
					converted_novel_pos_l.append(chr_s +"_"+ chr_pos_d[chr_s][j][0] +"_"+ chr_pos_d[chr_s][j][1])
#	print("converted_novel_pos_l:" ,converted_novel_pos_l, sep="\t")
	
	converted_info = ";".join(info_l[0:3]) +";"+ ",".join(converted_novel_pos_l)
#	print("converted_info: ", converted_info, sep="\t")
	
	if converted_info in novel_pos_dict:
		for i in range(len(line_l[1:])):
#			print(novel_pos_dict[converted_info][i], line_l[i+1])
			novel_pos_dict[converted_info][i] = str(int(novel_pos_dict[converted_info][i]) + int(line_l[i+1]))
#			print(novel_pos_dict[converted_info][i], line_l[i+1])
#		print("novel_pos_dict[converted_info]: ", novel_pos_dict[converted_info], sep="\t")
	else:
		novel_pos_dict[converted_info] = line_l[1:]
#		print("novel_pos_dict[converted_info]: ", novel_pos_dict[converted_info], sep="\t")
#	print()

for k,v in novel_pos_dict.items():
	print(k,"\t".join(v), sep="\t")

