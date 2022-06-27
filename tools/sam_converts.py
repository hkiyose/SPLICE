import sys
import re

def sam_convert(f):
	def get_start_end_on_genome(read_tmp):
		read_tmp_l = re.split('\t', read_tmp)
		cigar_str = read_tmp_l[5]
		cigar_l = cigar_str.split('N')

		map_starts = []
		map_ends = []
		i_count = 0

		for i in range(len(cigar_l)):
			i_count += 1#first exon
			cigar_type = [j for j in re.split(r'\d+', cigar_l[i]) if j != ""]
			cigar_len = [j for j in re.split(r'[a-zA-Z]+', cigar_l[i]) if j != ""]

			if i_count ==  1:
				map_start_pos = int(read_tmp_l[3])
				map_end = 0
				for i in range(0, len(cigar_type)):
	#				print(cigar_type[i], cigar_len[i])
					if cigar_type[i] == "M" or cigar_type[i] == "D":
						map_end += int(cigar_len[i])

				map_end_pos = map_start_pos + map_end -1
				map_starts.append(map_start_pos)
				map_ends.append(map_end_pos)

				if len(cigar_type) != len(cigar_len):#N cigar_len exists=>next exon exists
					nextExon_start_pos = map_end_pos + int(cigar_len[-1]) + 1
	#			print(nextExon_start_pos)
			else:#read has multiple exon
				map_start_pos = nextExon_start_pos

				map_end = 0
				for i in range(0, len(cigar_type)):
	#				print(cigar_type[i], cigar_len[i])
					if cigar_type[i] == "M" or cigar_type[i] == "D":
						map_end += int(cigar_len[i])

				map_end_pos = map_start_pos + map_end - 1
				map_starts.append(map_start_pos)
				map_ends.append(map_end_pos)

				if len(cigar_type) != len(cigar_len):
					nextExon_start_pos = map_end_pos + int(cigar_len[-1]) + 1
	#				print(nextExon_start_pos)
		return(map_starts,map_ends)

	def get_start_end_on_read(read_tmp):

		cigar_str = ""
		read_tmp_l = re.split('\t', read_tmp)
		cigar_type1 = re.split(r'\d+', read_tmp_l[5])
		cigar_type1 = re.split(r'\d+', read_tmp_l[5])
		cigar_len1 = re.split(r'[a-zA-Z]+', read_tmp_l[5])
			
		del cigar_type1[0]
		del cigar_len1[-1]

	#       print(cigar_type1)
	#       print(cigar_len1)
		
		if int(read_tmp_l[1]) & 0x10:#Unify strand
			cigar_type1.reverse()
			cigar_len1.reverse()
	#               print(cigar_type1)
	#               print(cigar_len1)
			
		for i in range(0, len(cigar_type1)):
			cigar_str += cigar_len1[i]+cigar_type1[i]
	#       print(cigar_str)
			
		cigar_l = cigar_str.split('N')
			
		read_starts = []
		read_ends = []
		read_start_pos = 0
		read_end_pos = 0
		i_count = 0

			
		for i in range(len(cigar_l)):
			i_count += 1#first exon
			cigar_type = [j for j in re.split(r'\d+', cigar_l[i]) if j != ""]
			cigar_len = [j for j in re.split(r'[a-zA-Z]+', cigar_l[i]) if j != ""]
					
			if i_count ==  1:
				for i in range(0, len(cigar_type)):
					if (cigar_type[i] == "S" or cigar_type[i] == "H") and i == 0:
						if cigar_type[i+1] == "I" or cigar_type[i+1] == "M":
							read_start_pos = int(cigar_len[i]) + 1
											
						else:#SやHの直後にDがある場合
							read_start_pos = int(cigar_len[i]) + int(cigar_len[i+1]) + 1
					elif cigar_type[i] == "I" or cigar_type[i] == "M":
						read_end_pos += int(cigar_len[i])
							
				read_end_pos += read_start_pos - 1
				read_starts.append(read_start_pos)
				read_ends.append(read_end_pos)
					
			else:#readが複数のexonを持つ場合
				for i in range(0, len(cigar_type)):
					if i == 0:
						if cigar_type[i] == "M" or cigar_type[i] == "I":
							read_start_pos = read_end_pos + 1
							read_end_pos += int(cigar_len[i])
											
						else :
							read_start_pos = read_end_pos + int(cigar_len[i]) + 1
									
					elif cigar_type[i] == "I" or cigar_type[i] == "M":
						read_end_pos += int(cigar_len[i])

				read_starts.append(read_start_pos)
				read_ends.append(read_end_pos)
			
		return(read_starts,read_ends)

	def convert_list_into_str(ls):
		str_l = map(str, ls)
		return ",".join(str_l)

	f1 = open(f)
	cs = ""
	read_id = ""
	for line in f1:
		if line[0] == "@":
			continue

		line = line.replace("\n", "")
		line_l = line.split("\t")
			
		if line_l[2] == "*":
			continue

	#	print(line_l)

		map_starts,map_ends = get_start_end_on_genome(line)
	#	print(map_starts,map_ends)
		read_starts,read_ends = get_start_end_on_read(line)
	#	print(read_starts,read_ends)

		if read_id == "":
			read_id = line_l[0]
			read_len = len(line_l[9])
		elif read_id != line_l[0]:
			read_id = line_l[0]
			read_len = len(line_l[9])
		elif read_id == line_l[0] and len(line_l[9]) != 1:
			read_len = len(line_l[9])
	#	print("read_len: ", read_len)

		map_starts_str = convert_list_into_str(map_starts)
		map_ends_str = convert_list_into_str(map_ends)
		read_starts_str = convert_list_into_str(read_starts)
		read_ends_str = convert_list_into_str(read_ends)
		
		for i in range(len(line_l)):
			if line_l[i].startswith("cs:"):
				cs = line_l[i]
		
		strand = "+"
		if int(line_l[1]) & 0x10:
			strand = "-"
			
		if  "/right" in line or "/left" in line:
			print(line_l[0], read_len, line_l[1], line_l[2], read_starts_str, read_ends_str, map_starts_str, map_ends_str, strand, line_l[4], cs, line_l[9], sep="\t")
								
		else:         
			print(line_l[0] + "/" + "None", read_len, line_l[1], line_l[2], read_starts_str, read_ends_str, map_starts_str, map_ends_str, strand, line_l[4], cs, line_l[9], sep="\t")
       

def sam_convert2(f, f_2, ol_rate):

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

				for k in range(len(cigar_s)):
					if cigar_s[k] == "M":
						match_count += 1
		return match_count

	f2 = open(f_2)#FASTQ
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
		elif line_num % 4 == 2:
			read_len_dict[read_name] = len(line)

	pre_read_id = ""
	pre_read_info = ""
	f1 = open(f) 
	for line in f1:
		line = line.replace("\n", "")
		line_l = line.split("\t")
		line = line.replace("\n", "")
		line_l_l = line_l[0].split("/")
		
		if pre_read_id != line_l_l[0] and pre_read_id != "":
			pre_read_info_l = pre_read_info.split(";")
			pre_read_info_l_l  = pre_read_info_l[0].split("/")
			print(pre_read_info_l_l[0],read_len_dict[pre_read_info_l_l[0]],pre_read_info_l[2],pre_read_info_l[3],pre_read_info_l[4],pre_read_info_l[5],pre_read_info_l[6],pre_read_info_l[7],pre_read_info_l[8],pre_read_info_l[9],pre_read_info_l[10],sep="\t") 
			pre_read_id = ""
			pre_read_info = ""
			
		if pre_read_id == "":
			pre_read_id = line_l_l[0]
			pre_read_info = line_l[0] +";"+ line_l[1] +";"+ line_l[2] +";"+ line_l[3] +";"+ line_l[4] +";"+ line_l[5] +";"+ line_l[6] +";"+ line_l[7] +";"+ line_l[8] +";"+ line_l[9] +";"+ line_l[10] 
			continue

		if pre_read_id == line_l_l[0]:
			
			if not "/left" in pre_read_info:
				if not "/right" in pre_read_info:
					pre_read_info_l = pre_read_info.split(";")
					pre_end_pos_l = re.split('[,/]', pre_read_info_l[5])
					softclip_right_add_len = int(pre_end_pos_l[-1])
				
			map_starts_l = re.split('[,/]', line_l[4])
			map_ends_l = re.split('[,/]', line_l[5])
			if len(line_l_l) > 3: 		
				if line_l_l[3] == "right":#Correct the position of the right softclip
					for i in range(len(map_starts_l)):
						map_starts_l[i] = int(map_starts_l[i]) + softclip_right_add_len
						map_ends_l[i] = int(map_ends_l[i]) + softclip_right_add_len
			map_starts_str = convert_list_into_str(map_starts_l)
			map_ends_str = convert_list_into_str(map_ends_l) 	
			
			read_region_l = [x for x in range(int(map_starts_l[0]), int(map_ends_l[-1])+1)] 
			
			pre_line_l = pre_read_info.split(";")
			pre_map_starts_l = re.split('[,/]', pre_line_l[4])
			pre_map_ends_l = re.split('[,/]', pre_line_l[5])
			pre_region_l = [x for x in range(int(pre_map_starts_l[0]), int(pre_map_ends_l[-1])+1)]
			
			read_region_set = set(read_region_l)
			pre_region_set = set(pre_region_l)
			read_match_count = match_count(line_l[10])
			pre_match_count = match_count(pre_line_l[10])
			
			and_readRegionSet_preRegionSet = read_region_set & pre_region_set
			
			overlap_eva = "None"
			overlap_rate = 0
			overlap_rate_pre = 0
			if len(and_readRegionSet_preRegionSet) > 0:#Avoid zero arithmetic
				overlap_rate = round(len(and_readRegionSet_preRegionSet)/len(read_region_set), 3)
				overlap_rate_pre = round(len(and_readRegionSet_preRegionSet)/len(pre_region_set), 3)
				overlap_eva = "Overlap"#If even one base is covered=>Overlap

			if overlap_rate <= float(ol_rate) and overlap_rate_pre <= float(ol_rate) or overlap_eva == "None":#Low overlap rate=>combine
				if int(read_region_l[0]) < int(pre_region_l[0]):#to the left
					pre_read_info = line_l[0] +";"+ line_l[1] +"/"+ pre_line_l[1]  +";"+ line_l[2] +"/"+ pre_line_l[2] +";"+ line_l[3] +"/"+ pre_line_l[3] +";"+ line_l[4] +"/"+ pre_line_l[4] +";"+ line_l[5] +"/"+ pre_line_l[5] +";"+ line_l[6] +"/"+ pre_line_l[6] +";"+ line_l[7] +"/"+ pre_line_l[7] +";"+ line_l[8] +"/"+ pre_line_l[8] +";"+ line_l[9] +"/"+ pre_line_l[9] +";"+ line_l[10] +"/"+ pre_line_l[10]
		
				elif int(read_region_l[0]) > int(pre_region_l[0]):#to the right
					pre_read_info = line_l[0] +";"+ pre_line_l[1] +"/"+ line_l[1] +";"+ pre_line_l[2] +"/"+ line_l[2] +";"+ pre_line_l[3] +"/"+ line_l[3] +";"+ pre_line_l[4] +"/"+ map_starts_str +";"+ pre_line_l[5] +"/"+ map_ends_str +";"+ pre_line_l[6] +"/"+ line_l[6] +";"+ pre_line_l[7] +"/"+ line_l[7] +";"+ pre_line_l[8] +"/"+ line_l[8] +";"+ pre_line_l[9] +"/"+ line_l[9] +";"+ pre_line_l[10] +"/"+ line_l[10]

			else:#overlap率が高い
				if read_match_count > pre_match_count:
					pre_read_info = line_l[0] +";"+ line_l[1] +";"+ line_l[2] +";"+ line_l[3] +";"+ line_l[4] +";"+ line_l[5] +";"+ line_l[6] +";"+ line_l[7] +";"+ line_l[8] +";"+ line_l[9] +";"+ line_l[10]
				elif read_match_count == pre_match_count:#If match_count is the same, give priority to the one with higher MAPQ.
					if line_l[9] > pre_line_l[9]:
						pre_read_info = line_l[0] +";"+ line_l[1] +";"+ line_l[2] +";"+ line_l[3] +";"+ line_l[4] +";"+ line_l[5] +";"+ line_l[6] +";"+ line_l[7] +";"+ line_l[8] +";"+ line_l[9] +";"+ line_l[10]

	pre_read_info_l = pre_read_info.split(";")
	pre_read_info_l_l  = pre_read_info_l[0].split("/")
	print(pre_read_info_l_l[0],read_len_dict[pre_read_info_l_l[0]],pre_read_info_l[2],pre_read_info_l[3],pre_read_info_l[4],pre_read_info_l[5],pre_read_info_l[6],pre_read_info_l[7],pre_read_info_l[8],pre_read_info_l[9],pre_read_info_l[10],sep="\t") 

def sam_convert3(f, f_2, min_sc_len):

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

				for k in range(len(cigar_s)):
					if cigar_s[k] == "M":
						match_count += 1
		return match_count

	f2 = open(f_2)#FASTQ
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
		elif line_num % 4 == 2:
			read_len_dict[read_name] = len(line)

	pre_read_id = ""
	pre_read_info = ""
	f1 = open(f) 
	for line in f1:
		line = line.replace("\n", "")
		line_l = line.split("\t")
		line = line.replace("\n", "")
		line_l_l = line_l[0].split("/")
		
		if pre_read_id != line_l_l[0] and pre_read_id != "":#output
			pre_read_info_l = pre_read_info.split(";")
			pre_read_info_l_l  = pre_read_info_l[0].split("/")
			print(pre_read_info_l_l[0],read_len_dict[pre_read_info_l_l[0]],pre_read_info_l[2],pre_read_info_l[3],pre_read_info_l[4],pre_read_info_l[5],pre_read_info_l[6],pre_read_info_l[7],pre_read_info_l[8],pre_read_info_l[9],pre_read_info_l[10],sep="\t") 
			pre_read_id = ""
			pre_read_info = ""
			
		if pre_read_id == "":
			pre_read_id = line_l_l[0]
			pre_read_info = line_l[0] +";"+ line_l[1] +";"+ line_l[2] +";"+ line_l[3] +";"+ line_l[4] +";"+ line_l[5] +";"+ line_l[6] +";"+ line_l[7] +";"+ line_l[8] +";"+ line_l[9] +";"+ line_l[10] 
			continue

		if pre_read_id == line_l_l[0]:
			
			if not "/left" in pre_read_info:
				if not "/right" in pre_read_info:
					pre_read_info_l = pre_read_info.split(";")
					pre_end_pos_l = re.split('[,/]', pre_read_info_l[5])
					softclip_right_add_len = int(pre_end_pos_l[-1])
				
			map_starts_l = re.split('[,/]', line_l[4])
			map_ends_l = re.split('[,/]', line_l[5])
			if len(line_l_l) == 2: 		
				if line_l_l[1] == "right":
					for i in range(len(map_starts_l)):
	#					print(line_l)
						map_starts_l[i] = int(map_starts_l[i]) + softclip_right_add_len + 1
						map_ends_l[i] = int(map_ends_l[i]) + softclip_right_add_len + 1
			map_starts_str = convert_list_into_str(map_starts_l)
			map_ends_str = convert_list_into_str(map_ends_l) 	
			
			read_region_l = [x for x in range(int(map_starts_l[0]), int(map_ends_l[-1])+1)] 
			
			pre_line_l = pre_read_info.split(";")
			pre_map_starts_l = re.split('[,/]', pre_line_l[4])
			pre_map_ends_l = re.split('[,/]', pre_line_l[5])
			pre_region_l = [x for x in range(int(pre_map_starts_l[0]), int(pre_map_ends_l[-1])+1)]
			read_region_set = set(read_region_l)
			pre_region_set = set(pre_region_l)
			read_match_count = match_count(line_l[10])
			pre_match_count = match_count(pre_line_l[10])
			
			and_readRegionSet_preRegionSet = read_region_set & pre_region_set
			
			overlap_eva = "None"
			overlap_rate = 0
			overlap_rate_pre = 0
			if len(and_readRegionSet_preRegionSet) > 0:
				overlap_rate = round(len(and_readRegionSet_preRegionSet)/len(read_region_set), 3)
				overlap_rate_pre = round(len(and_readRegionSet_preRegionSet)/len(pre_region_set), 3)
				overlap_eva = "exist"

			if overlap_rate <= float(min_sc_len) and overlap_rate_pre <= float(min_sc_len) or overlap_eva == "None":
				if int(read_region_l[0]) < int(pre_region_l[0]):
					pre_read_info = line_l[0] +";"+ line_l[1] +"/"+ pre_line_l[1]  +";"+ line_l[2] +"/"+ pre_line_l[2] +";"+ line_l[3] +"/"+ pre_line_l[3] +";"+ line_l[4] +"/"+ pre_line_l[4] +";"+ line_l[5] +"/"+ pre_line_l[5] +";"+ line_l[6] +"/"+ pre_line_l[6] +";"+ line_l[7] +"/"+ pre_line_l[7] +";"+ line_l[8] +"/"+ pre_line_l[8] +";"+ line_l[9] +"/"+ pre_line_l[9] +";"+ line_l[10] +"/"+ pre_line_l[10]
		
				elif int(read_region_l[0]) > int(pre_region_l[0]):
					pre_read_info = line_l[0] +";"+ pre_line_l[1] +"/"+ line_l[1] +";"+ pre_line_l[2] +"/"+ line_l[2] +";"+ pre_line_l[3] +"/"+ line_l[3] +";"+ pre_line_l[4] +"/"+ map_starts_str +";"+ pre_line_l[5] +"/"+ map_ends_str +";"+ pre_line_l[6] +"/"+ line_l[6] +";"+ pre_line_l[7] +"/"+ line_l[7] +";"+ pre_line_l[8] +"/"+ line_l[8] +";"+ pre_line_l[9] +"/"+ line_l[9] +";"+ pre_line_l[10] +"/"+ line_l[10]

			else:
				if line_l[9] > pre_line_l[9]:
					pre_read_info = line_l[0] +";"+ line_l[1] +";"+ line_l[2] +";"+ line_l[3] +";"+ line_l[4] +";"+ line_l[5] +";"+ line_l[6] +";"+ line_l[7] +";"+ line_l[8] +";"+ line_l[9] +";"+ line_l[10]
				elif line_l[9] == pre_line_l[9]:
					if read_match_count > pre_match_count:
						pre_read_info = line_l[0] +";"+ line_l[1] +";"+ line_l[2] +";"+ line_l[3] +";"+ line_l[4] +";"+ line_l[5] +";"+ line_l[6] +";"+ line_l[7] +";"+ line_l[8] +";"+ line_l[9] +";"+ line_l[10]

	pre_read_info_l = pre_read_info.split(";")
	pre_read_info_l_l  = pre_read_info_l[0].split("/")
	print(pre_read_info_l_l[0],read_len_dict[pre_read_info_l_l[0]],pre_read_info_l[2],pre_read_info_l[3],pre_read_info_l[4],pre_read_info_l[5],pre_read_info_l[6],pre_read_info_l[7],pre_read_info_l[8],pre_read_info_l[9],pre_read_info_l[10],sep="\t") 

def sam_convert4(f):

	f1 = open(f)
	for line in f1:
		line = line.replace("\n", "")
		line_l = line.split("\t")
#		print(line)
		if "/" in line_l[2]:
			continue	
		for i in range(len(line_l)):
			if line_l[i].startswith("cs:"):
				cs = line_l[i]
#		print(cs)
	
		cs_l = cs.split('~')
		cigar_s = ""
			
		for i in range(len(cs_l)):
			cs_l_type = re.findall('[:*+-]', cs_l[i])
			cs_l_len = re.split('[:*+-]', cs_l[i])	
			if i == 0:
				del cs_l_type[0:2]
				del cs_l_len[0:3]
			else:
				del cs_l_len[0]
			
#			print("cs_l_type: ", len(cs_l_type), cs_l_type )
#			print("cs_l_len: ", len(cs_l_len), cs_l_len)
	
			for j in range(len(cs_l_type)):
				if cs_l_type[j] == ":":
					cigar_s += "M" * int(cs_l_len[j])
				elif cs_l_type[j] == "*":
					cigar_s += "R"
				elif cs_l_type[j] == "+":
					cigar_s += "I" * len(cs_l_len[j])
				else:#cs_l_type[j] == "-"
					cigar_s += "D" * len(cs_l_len[j])
#		print(cigar_s)	
		print(cigar_s.count("M"),cigar_s.count("M")/int(line_l[1]), line, cigar_s, sep="\t")	
#		print("")
