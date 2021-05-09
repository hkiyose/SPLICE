import sys

f = open(sys.argv[1])
num_line = 0
dict_qscore = {}

for line in f:
	num_line += 1
	line = line.replace("\n", "")
	sum_q = 0
	ave_q = 0
#	print(line)
	if num_line%4 == 1:
		line_l = line.split()
		read_name = line
	elif num_line%4 == 2:
		seq = line
	elif num_line%4 == 0:
		for q in line:
			q2 = ord(q) - 33
			sum_q += q2
		ave_q = int(sum_q/len(line))
		if ave_q >= int(sys.argv[2]):
			print(read_name)
			print(seq)
			print("+")
			print(line)
