import sys
import re

base_convert_dict = {"A":"T", "T":"A", "G":"C", "C":"G"}
def get_rev_seq(seq):
	seq2 = ""
	for i in range(len(seq)):
		seq2 += base_convert_dict[seq[i]]
	seq2 = seq2[::-1]
	return seq2
	

f = open(sys.argv[1])#refseq_exonNum
gene_strand_dict = {}
for line in f:
	line = line.replace("\n", "")
	line_l = line.split("\t")
#	print(line_l[12], line_l[3], sep="\t")
	if not line_l[12] in gene_strand_dict:
		gene_strand_dict[line_l[12]] = line_l[3]

#for k,v in gene_strand_dict.items():
#	print(k, v, sep="\t")

f2 = open(sys.argv[2])#
read_3prime_seq_dict = {}#read_3prime_seq_dict[read_id] = start_pos, end_posi, strand...Invert the array when "-".
for line in f2:
	line = line.replace("\n", "")
	line_l = line.split("\t")
	if "/" in line_l[11] or "novel" in line_l[11]:
		continue
#	print(line)
#	print("line_l[11]: ", line_l[11], sep="\t")
#	print("gene_strand_dict[line_l[11]]: ", gene_strand_dict[line_l[11]], sep="\t")

	start_pos_l = line_l[4].split(",")
	end_pos_l = line_l[5].split(",")	
	
	if gene_strand_dict[line_l[11]] == "+":
		v = [str(int(end_pos_l[-1])+1), line_l[1], line_l[8]]
		read_3prime_seq_dict[line_l[0]] = v
		
	else:
		v = ["0", str(int(start_pos_l[0])-1), line_l[8]]
		read_3prime_seq_dict[line_l[0]] = v
	
#	print("read_3prime_seq_dict[line_l[0]]: ", read_3prime_seq_dict[line_l[0]], sep="\t")
	#['327', '362', '-']	
#	print("")

f3 = open(sys.argv[3])#fastq
count = 0
pA_count_dict = {}
pA_read_dict = {}
for line in f3:
	line = line.replace("\n", "")
	count += 1
#	print(str(count % 4), line, sep="\t")
	if count % 4 == 1:
		line_l = line.split()
		read_id = line_l[0][1:]
##		print("read_id: ", read_id, sep="\t")
	
	elif count % 4 == 2:
		if read_id in read_3prime_seq_dict:
##			print("read_3prime_seq_dict[read_id]: ", read_3prime_seq_dict[read_id], sep="\t")
			sc_3prime_seq = ""
			if read_3prime_seq_dict[read_id][2] == "+":
#				print(line)
##				print(line[int(read_3prime_seq_dict[read_id][0]):int(read_3prime_seq_dict[read_id][1])])
				sc_3prime_seq = line[int(read_3prime_seq_dict[read_id][0]):int(read_3prime_seq_dict[read_id][1])]
			else:
#				print("rev: ", get_rev_seq(line))
				rev_seq = get_rev_seq(line)
##				print(rev_seq[int(read_3prime_seq_dict[read_id][0]):int(read_3prime_seq_dict[read_id][1])])
				sc_3prime_seq = rev_seq[int(read_3prime_seq_dict[read_id][0]):int(read_3prime_seq_dict[read_id][1])]

			pA_eva = False
##			print(int(int(sys.argv[4]) * float(sys.argv[5])))
			for i in range(len(sc_3prime_seq)-int(sys.argv[4])+1):
##				print(sc_3prime_seq[i:i+int(sys.argv[4])])
				seg_seq = sc_3prime_seq[i:i+int(sys.argv[4])]	
				if read_3prime_seq_dict[read_id][0] == "0":
					if seg_seq.count("T") >= int(int(int(sys.argv[4]) * float(sys.argv[5]))):
##						print("polyT: ", seg_seq.count("T"))	
						pA_eva = True
				else:		
					if seg_seq.count("A") >= int(int(int(sys.argv[4]) * float(sys.argv[5]))):
##						print("polyA: ", seg_seq.count("A"))
						pA_eva = True
##			print("")
			if pA_eva == True:
#				print(sc_3prime_seq)			
				pA_read_dict[read_id] = 0
#for k,v in sorted(pA_count_dict.items(), key=lambda x:x[0]):
#	print(k,v,sep="\t")			
		
#for k,v in pA_read_dict.items():
#	print(k,v,sep="\t")

f4 = open(sys.argv[6])
for line in f4:
	line = line.replace("\n","")
	line_l = line.split("\t")
#	print(line)
	if "KNOWN" in line_l[0]:
		continue	
	read_id_l = line_l[-1].split(",")
	read_id_l2 = []
	for i in range(len(read_id_l)):
		if read_id_l[i] in pA_read_dict:
			read_id_l2.append(read_id_l[i])
	if len(read_id_l2) >= 1:
		print(line_l[0], len(read_id_l2), line_l[2], line_l[3], line_l[4], line_l[5], line_l[6], line_l[7], line_l[8], ",".join(read_id_l2), sep="\t")
#		print(line_l[0], len(read_id_l2), "*", line_l[3], line_l[4], line_l[5], line_l[6], line_l[7], line_l[8], line_l[9], ",".join(read_id_l2), sep="\t")

