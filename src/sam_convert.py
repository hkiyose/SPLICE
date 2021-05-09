import sys
import re

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

	
f = open(sys.argv[1])
	
cs = ""
read_id = ""
for line in f:
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
       

