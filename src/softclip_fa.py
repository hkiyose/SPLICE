import sys
import re 

f = open(sys.argv[1])
map_pos_on_read_dict = {}
for line in f:
	line = line.replace("\n", "")
	line_l = line.split("\t")
#	print(line_l)
	map_start_pos_on_read_l = re.split('[,/]', line_l[4])
	map_end_pos_on_read_l = re.split('[,/]', line_l[5])
	start_pos = map_start_pos_on_read_l[0]
	end_pos =  map_end_pos_on_read_l[-1]
	map_pos_on_read_dict[line_l[0]] = start_pos, end_pos

#for k,v in map_pos_on_read_dict.items():
#	print(k, v)

f2 = open(sys.argv[2])#FASTQ_181213
line_num = 0
for line in f2:
	line_num += 1
	line = line.replace("\n", "")
	line_l = line.split()
#	print(line_num, line)
	
	if line_num % 4 == 1:
		read_id = line_l[0]
		read_id = read_id[1:]
			
		start_pos = ""
		end_pos = ""
		if read_id in map_pos_on_read_dict:
			start_pos, end_pos = map_pos_on_read_dict[read_id]
#			print(line_l[0], start_pos, end_pos)

	if line_num % 4 == 2:
		seq = line
		if not start_pos == "":
#			left_seq = seq[0:int(start_pos)-1]
			left_seq = seq[0:int(start_pos)]
			left_seq = left_seq[:-1]
			right_seq = seq[int(end_pos)+1:len(seq)]
#			print("left_seq: ", left_seq, sep="\t")
#			print("right_seq: ", right_seq, sep="\t")
			if len(left_seq) >= int(sys.argv[3]):
				print(">" + read_id + "/left")	
				print(left_seq)
			if len(right_seq) >= int(sys.argv[3]):
				print(">" + read_id + "/right")
				print(right_seq)
#		print("")

