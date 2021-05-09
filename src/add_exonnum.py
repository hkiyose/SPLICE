import sys

def remove_duplicates(x):
	y=[]
#	dup_count = 0
	for i in x:
		if i not in y:
			y.append(i)
	return y

def get_uniq_sort_dict(dict):
	for k,v in dict.items():
#		print(k,v,sep="\t")
#		print(len(v))
		uniq_v = remove_duplicates(v)
#		print("uniq_v: ", uniq_v, sep="\t")
		uniq_v_sort = sorted(uniq_v)
#		print("uniq_v_sort: ", uniq_v_sort, sep="\t")
		exon_dict_uniq[k] = uniq_v_sort
#		print("exon_dict_uniq[k]: ", len(exon_dict_uniq[k]), exon_dict_uniq[k], sep="\t")
	return exon_dict_uniq

def get_exon_num(dict):
#	print(dict)
	for k,v in dict.items():
		exon_num = {}
#		print(k,v,sep="\t")
		for i in range(len(v)):
#			print(k, v[i][0], v[i][1])
			exon_pos = v[i][0],v[i][1]
			exon_num[exon_pos] = i+1
			exon_dict_uniq_num[k] = exon_num
	return exon_dict_uniq_num

f1 = open(sys.argv[1])
f2 = open(sys.argv[1])

gene_dict = {}
exon_dict = {}
exon_dict_uniq = {}
exon_dict_uniq_num = {}
pre_gene_name = ""
gene_list = []
uniq_dict = {}

for line in f1:
	line = line.replace("\n","")
	line_l = line.split("\t")
	if line.startswith("#"):
		continue
	exon_start = line_l[9].split(",")
	exon_end = line_l[10].split(",")
	del exon_start[-1]
	del exon_end[-1]

	for i in range(len(exon_start)):
		exon_dict.setdefault(line_l[12], []).append([int(exon_start[i]),int(exon_end[i])])
#print(exon_dict)
exon_uniq_sort_dict = get_uniq_sort_dict(exon_dict)
#print(exon_uniq_sort_dict)
exon_uniq_sort_num_dict = get_exon_num(exon_uniq_sort_dict)
#print(exon_uniq_sort_num_dict)

for line in f2:
	line = line.replace("\n","")
	line_l = line.split("\t")
	if line.startswith("#"):
		continue
	exon_start = line_l[9].split(",")
	exon_end = line_l[10].split(",")
	del exon_start[-1]
	del exon_end[-1]

	exon_num_list = []
	for i in range(len(exon_start)):
		exon_pos = int(exon_start[i]),int(exon_end[i])
		exon_num_list.append(exon_pos)
#    print(exon_num_list)
	exon_num_iso_list = []
	exon_num_iso_str = ""
	for exon in exon_num_list:
		exon_num_iso_list.append(int(exon_uniq_sort_num_dict[str(line_l[12])][exon]))
#    print(exon_num_iso_list)
	mapped_list = map(str, exon_num_iso_list)
	mapped_list_str = ",".join(mapped_list)
#    print(mapped_list_str)
	print(line+"\t"+mapped_list_str)
