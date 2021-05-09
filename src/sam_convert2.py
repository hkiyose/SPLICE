import sys
import re

def convert_list_into_str(ls):
        str_l = map(str, ls)
        return ",".join(str_l)

def match_count(cs):
	match_count = 0
	cs_l = cs.split("/")
	for i in range(len(cs_l)):
		cs_l_l = cs_l[i].split('~')

		for j in range(len(cs_l_l)):
			cs_l_type = re.findall('[:*+-]', cs_l_l[j])
			cs_l_len = re.split('[:*+-]', cs_l_l[j])

			if j == 0:
				del cs_l_type[0:2]
				del cs_l_len[0:3]
			else:
				del cs_l_len[0]
		
#			print("cs_l_type: ",len(cs_l_type), cs_l_type)
#			print("cs_l_len: ", len(cs_l_len), cs_l_len)

			cigar_s = ""		                        
			for k in range(len(cs_l_type)):
				if cs_l_type[k] == ":":
					cigar_s += "M" * int(cs_l_len[k])
				elif cs_l_type[k] == "*":
					cigar_s += "R"
				elif cs_l_type[k] == "+":
					cigar_s += "I" * len(cs_l_len[k])
				else:#cs_l_type[k] == "-"
					cigar_s += "D" * len(cs_l_len[k])

#			print("cigar_s: ", cigar_s)
			for k in range(len(cigar_s)):
				if cigar_s[k] == "M":
					match_count += 1
#			print("match_count: ", match_count)
	return match_count

f2 = open(sys.argv[2])#FASTQ
read_name = ""
read_len_dict = {}
line_num = 0
for line in f2:
	line_num += 1
	line = line.replace("\n", "")
	line_l = line.split()
	line_l_l = line_l[0].split("/")
	if line_num % 4 == 1:
		read_name = line_l_l[0]
		read_name = read_name[1:]	
#		print(read_name)	
	elif line_num % 4 == 2:
		read_len_dict[read_name] = len(line)
#for k,v in read_len_dict.items():
#	print(k, v)	

pre_read_id = ""
pre_read_info = ""
f = open(sys.argv[1]) 
for line in f:
	line = line.replace("\n", "")
	line_l = line.split("\t")
##	print(">>", line)
	line = line.replace("\n", "")
	line_l_l = line_l[0].split("/")
	
	if pre_read_id != line_l_l[0] and pre_read_id != "":
#		print(pre_read_info)
		pre_read_info_l = pre_read_info.split(";")
		pre_read_info_l_l  = pre_read_info_l[0].split("/")
		print(pre_read_info_l_l[0],read_len_dict[pre_read_info_l_l[0]],pre_read_info_l[2],pre_read_info_l[3],pre_read_info_l[4],pre_read_info_l[5],pre_read_info_l[6],pre_read_info_l[7],pre_read_info_l[8],pre_read_info_l[9],pre_read_info_l[10],sep="\t") 
		pre_read_id = ""
		pre_read_info = ""
##		print()
		
	if pre_read_id == "":
		pre_read_id = line_l_l[0]
		pre_read_info = line_l[0] +";"+ line_l[1] +";"+ line_l[2] +";"+ line_l[3] +";"+ line_l[4] +";"+ line_l[5] +";"+ line_l[6] +";"+ line_l[7] +";"+ line_l[8] +";"+ line_l[9] +";"+ line_l[10] 
		continue

	if pre_read_id == line_l_l[0]:
		
#		if pre_read_info == "":
#			pre_read_info = line_l[0] +":"+ line_l[1] +":"+ line_l[2] +":"+ line_l[3] +":"+ line_l[4] +":"+ line_l[5] +":"+ line_l[6] +":"+ line_l[7] +":"+ line_l[8] +":"+ line_l[9] +":"+ line_l[10] 
		
		if not "/left" in pre_read_info:
			if not "/right" in pre_read_info:
				pre_read_info_l = pre_read_info.split(";")
				pre_end_pos_l = re.split('[,/]', pre_read_info_l[5])
				softclip_right_add_len = int(pre_end_pos_l[-1])
			
		map_starts_l = re.split('[,/]', line_l[4])
		map_ends_l = re.split('[,/]', line_l[5])
#		print(map_starts_l)
#		print(map_ends_l)
		if len(line_l_l) > 3: 		
			if line_l_l[3] == "right":#Correct the position of the right softclip
				for i in range(len(map_starts_l)):
#					print(line_l)
					map_starts_l[i] = int(map_starts_l[i]) + softclip_right_add_len
					map_ends_l[i] = int(map_ends_l[i]) + softclip_right_add_len
		map_starts_str = convert_list_into_str(map_starts_l)
		map_ends_str = convert_list_into_str(map_ends_l) 	
##		print("map_starts_l: ", map_starts_l)
##		print("map_ends_l: ", map_ends_l)
		
		read_region_l = [x for x in range(int(map_starts_l[0]), int(map_ends_l[-1])+1)] 
##		print("read_region_l: ", read_region_l)
		
		pre_line_l = pre_read_info.split(";")
#		print(">", pre_line_l)
		pre_map_starts_l = re.split('[,/]', pre_line_l[4])
		pre_map_ends_l = re.split('[,/]', pre_line_l[5])
#		print(pre_line_l)
#		print("pre_map_starts_l: ", pre_map_starts_l)
#		print("pre_map_ends_l: ", pre_map_ends_l)
#		print("pre_map_starts_l[0]: ", pre_map_starts_l[0])
#		print("pre_map_ends_l[-1]: ", pre_map_ends_l[-1])
		pre_region_l = [x for x in range(int(pre_map_starts_l[0]), int(pre_map_ends_l[-1])+1)]
##		print("pre_region_l: ", pre_region_l)
		
		read_region_set = set(read_region_l)
		pre_region_set = set(pre_region_l)
		read_match_count = match_count(line_l[10])
		pre_match_count = match_count(pre_line_l[10])
##		print("read_match_count: ", read_match_count)
##		print("pre_match_count: ", pre_match_count)
#		print("read_region_set: ", read_region_set)
#		print("pre_region_set: ", pre_region_set)
		
		and_readRegionSet_preRegionSet = read_region_set & pre_region_set
		
		overlap_eva = "None"
		overlap_rate = 0
		overlap_rate_pre = 0
		if len(and_readRegionSet_preRegionSet) > 0:#Avoid zero arithmetic
			overlap_rate = round(len(and_readRegionSet_preRegionSet)/len(read_region_set), 3)
			overlap_rate_pre = round(len(and_readRegionSet_preRegionSet)/len(pre_region_set), 3)
			overlap_eva = "Overlap"#If even one base is covered=>Overlap

#		non_overlap_set = pre_region_set - read_region_set
#		print("non_overlap_set: ", non_overlap_set)
#		print("non_overlap_len: ", len(non_overlap_set))
		
##		print("overlap_rate: ", overlap_rate)
##		print("overlap_rate_pre: ", overlap_rate_pre)		
##		print("overlap_eva: ", overlap_eva)
		if overlap_rate <= float(sys.argv[3]) and overlap_rate_pre <= float(sys.argv[3]) or overlap_eva == "None":#Low overlap rate=>combine
#		if len(non_overlap_set) >= int(sys.argv[3]):
##			print("low overlap")
#			print(read_region_l)
#			print(pre_region_l)
#			print(line_l_l)
			if int(read_region_l[0]) < int(pre_region_l[0]):#to the left
				pre_read_info = line_l[0] +";"+ line_l[1] +"/"+ pre_line_l[1]  +";"+ line_l[2] +"/"+ pre_line_l[2] +";"+ line_l[3] +"/"+ pre_line_l[3] +";"+ line_l[4] +"/"+ pre_line_l[4] +";"+ line_l[5] +"/"+ pre_line_l[5] +";"+ line_l[6] +"/"+ pre_line_l[6] +";"+ line_l[7] +"/"+ pre_line_l[7] +";"+ line_l[8] +"/"+ pre_line_l[8] +";"+ line_l[9] +"/"+ pre_line_l[9] +";"+ line_l[10] +"/"+ pre_line_l[10]
	
			elif int(read_region_l[0]) > int(pre_region_l[0]):#to the right
				pre_read_info = line_l[0] +";"+ pre_line_l[1] +"/"+ line_l[1] +";"+ pre_line_l[2] +"/"+ line_l[2] +";"+ pre_line_l[3] +"/"+ line_l[3] +";"+ pre_line_l[4] +"/"+ map_starts_str +";"+ pre_line_l[5] +"/"+ map_ends_str +";"+ pre_line_l[6] +"/"+ line_l[6] +";"+ pre_line_l[7] +"/"+ line_l[7] +";"+ pre_line_l[8] +"/"+ line_l[8] +";"+ pre_line_l[9] +"/"+ line_l[9] +";"+ pre_line_l[10] +"/"+ line_l[10]

		else:#overlap率が高い
##			print("high overlap")#Give priority to the one with the highest match_count.
			if read_match_count > pre_match_count:
##				print("read_match_count > pre_match_count")
				pre_read_info = line_l[0] +";"+ line_l[1] +";"+ line_l[2] +";"+ line_l[3] +";"+ line_l[4] +";"+ line_l[5] +";"+ line_l[6] +";"+ line_l[7] +";"+ line_l[8] +";"+ line_l[9] +";"+ line_l[10]
			elif read_match_count == pre_match_count:#If match_count is the same, give priority to the one with higher MAPQ.
				if line_l[9] > pre_line_l[9]:
##					print("line_l[9] > pre_line_l[9]")
					pre_read_info = line_l[0] +";"+ line_l[1] +";"+ line_l[2] +";"+ line_l[3] +";"+ line_l[4] +";"+ line_l[5] +";"+ line_l[6] +";"+ line_l[7] +";"+ line_l[8] +";"+ line_l[9] +";"+ line_l[10]
				
				
##	print("")

pre_read_info_l = pre_read_info.split(";")
pre_read_info_l_l  = pre_read_info_l[0].split("/")
print(pre_read_info_l_l[0],read_len_dict[pre_read_info_l_l[0]],pre_read_info_l[2],pre_read_info_l[3],pre_read_info_l[4],pre_read_info_l[5],pre_read_info_l[6],pre_read_info_l[7],pre_read_info_l[8],pre_read_info_l[9],pre_read_info_l[10],sep="\t") 
