import sys
import re
from collections import OrderedDict

def annot_exonnum(f, f_2, mq_filt, min_gap_1, min_gap_2):
	f1 = open(f)

	exon_dict = OrderedDict()
	exon_dict_sort = OrderedDict()
	RNAid2name = OrderedDict()
	exon_pos2RNAid  = OrderedDict()
	exon_pos2RNAid_sort_dict = OrderedDict()

	for line in f1:

		line = line.replace("\n", "")
		line_l= line.split("\t")
	#	print(line_l[9])
		if not len(line_l) >= 10:
			continue 
		exon_start = line_l[9].split(",")
		exon_end = line_l[10].split(",")
		exon_num = line_l[16].split(",")

		del exon_start[-1]
		del exon_end[-1]

		RNAid2name[line_l[1]] = line_l[12]

		for i in range(len(exon_start)):
			exon_dict.setdefault(line_l[2], []).append([str(int(exon_start[i])+1), exon_end[i]])
			exon_position = line_l[2] + "\t" + str(int(exon_start[i])+1) + "\t" + exon_end[i]
			if not exon_position in exon_pos2RNAid:
				exon_pos2RNAid[exon_position] = []
				exon_pos2RNAid[exon_position].append([line_l[12],exon_num[i]])
			else:
				if not [line_l[12],exon_num[i]] in exon_pos2RNAid[exon_position]:
					exon_pos2RNAid[exon_position].append([line_l[12],exon_num[i]])

	#exon_dictの重複の削除とsort
	for chr in exon_dict:
		exon_dict_chr = exon_dict[chr]
		exon_dict_chr = list(map(list, set(map(tuple, exon_dict_chr))))
		exon_dict_chr.sort(key=lambda x:int(x[0]))
		exon_dict_chr.sort(key=lambda x:int(x[1]))
		exon_dict_sort[chr] = exon_dict_chr
	#print(exon_dict_sort)
	exon_pos2RNAid['novel\t*\t*'] = ('novel', '*')

	exon_pos_tmp_l = []
	exon_start_pos_tmp = ""
	seg_exon_dict = {}
	chr_seg_region_dict = {}
	for k,v in exon_dict_sort.items():
	#	print(k,len(v),v,sep="\t")
		chr_seg_region_dict[k] = []
		n = 0
		pos_l = v
		for i in range(len(pos_l)):
	##		print("pos_l[i]: ", pos_l[i], sep="\t")
			
			n += 1
			if n == 1:	
	#			print(pos_l[i][0])
				exon_start_pos_tmp = pos_l[i][0]
	#			print("exon_start_pos_tmp: ", exon_start_pos_tmp, sep="\t")
				exon_pos_tmp_l.append(pos_l[i])
	#			print("exon_pos_tmp_l: ", exon_pos_tmp_l, sep="\t")

	#		print(pos_l[i][0], pos_l[i][1], sep="\t")
			exon_pos_tmp_l.append(pos_l[i])
	##		print("exon_pos_tmp_l: ", exon_pos_tmp_l, sep="\t")

			if n == 100:
	##			print("len(exon_pos_tmp_l): ",len(exon_pos_tmp_l), sep="\t")
	#			print(pos_l[i][1])
				seg_k = (k, exon_start_pos_tmp ,pos_l[i][1])
				seg_k2 = (exon_start_pos_tmp ,pos_l[i][1])
	##			print("seg_k: ", seg_k, sep="\t")
	##			print("seg_v: ", exon_pos_tmp_l, sep="\t")
				seg_exon_dict[seg_k] = exon_pos_tmp_l	
				chr_seg_region_dict[k].append(list(seg_k2))
	##			print("seg_exon_dict[seg_k]: ", seg_exon_dict[seg_k], sep="\t")
	##			print("chr_seg_region_dict[k]: ", chr_seg_region_dict[k], sep="\t")
	##			print("seg_exon_dict: ", seg_exon_dict, sep="\t")	
	#			for k,v in seg_exon_dict.items():
	#				print(k,v,sep="\t")

				n = 0
				exon_pos_tmp_l = []
	##			print()

	##	print("len(exon_pos_tmp_l): ",len(exon_pos_tmp_l), sep="\t")	
		seg_k = (k, exon_start_pos_tmp ,pos_l[i][1])
		seg_k2 = (exon_start_pos_tmp ,pos_l[i][1])
	##	print("seg_k: ", seg_k, sep="\t")
	##	print("seg_v: ", exon_pos_tmp_l, sep="\t")
		seg_exon_dict[seg_k] = exon_pos_tmp_l
		chr_seg_region_dict[k].append(list(seg_k2))
	##	print("seg_exon_dict[seg_k]: ", seg_exon_dict[seg_k], sep="\t")
	##	print("chr_seg_region_dict[k]: ", chr_seg_region_dict[k], sep="\t")
		n = 0
		exon_pos_tmp_l = []


	##	print()
		

	#print(seg_exon_dict)
	#for k,v in seg_exon_dict:
	#	for i in v:
	#		print(k, v[i], sep="\t")
	#	print(k, len(v), v, sep="\t")

	segment_f = open(f_2)
	for line in segment_f:
		line = line.replace("\n", "")
		line_l = line.split("\t")
	##	print(">>>", line)
		annotated_exon_position = []
		uniq_exon_position = []
		annot = []
		read_condition = ""
		read_conditions = ""
		annot_gene_exonNum = ""
		annot_gene_name = ""
				
	#	if "chrM" in line_l[3]:
	#		continue
		
		map_starts_l = line_l[6].split("/")
		map_ends_l = line_l[7].split("/") 
		mq_l = line_l[9].split("/")
		chr_l = line_l[3].split("/")

		i_count = 0	
		for i in range(len(map_starts_l)):
			i_count += 1
	#		print("chr_l[i]: ", chr_l[i])
			
			if int(mq_l[i][0]) < int(mq_filt):#mq_filter
				continue

			if not chr_l[i] in exon_dict_sort:
				annot = ["ND"]
				gene_annot = ["ND"]

			else:
				novel_count = 0

				map_starts_l_l = map_starts_l[i].split(",")
				map_ends_l_l = map_ends_l[i].split(",") 

				j_count = 0
	##			print("chr_l[i]: ", chr_l[i], sep="\t")
				for j in range(len(map_starts_l_l)):
					j_count += 1
	##				print("map_starts_l_l[j]: ", map_starts_l_l[j])
	##				print("map_ends_l_l[j]: ", map_ends_l_l[j])

					highest_overlap_score = -1000000#
	##				print(chr_seg_region_dict)	
					for l in range(len(chr_seg_region_dict[chr_l[i]])):
	##					print(chr_seg_region_dict[chr_l[i]][l])
						if int(chr_seg_region_dict[chr_l[i]][l][0]) <= int(map_starts_l_l[j]) <= int(chr_seg_region_dict[chr_l[i]][l][1]) or int(chr_seg_region_dict[chr_l[i]][l][0]) <= int(map_ends_l_l[j]) <= int(chr_seg_region_dict[chr_l[i]][l][1]):		
	##						print("MATCH!!!")
	#						for exon in exon_dict_sort[chr_l[i]]:
							for exon in seg_exon_dict[(chr_l[i], chr_seg_region_dict[chr_l[i]][l][0], chr_seg_region_dict[chr_l[i]][l][1])]:
								if int(exon[0]) <= int(map_starts_l_l[j]) <= int(exon[1]) or int(exon[0]) <= int(map_ends_l_l[j]) <= int(exon[1]) or (int(map_starts_l_l[j]) < int(exon[0]) and int(exon[1]) < int(map_ends_l_l[j])):
	##								print("ref_start: ", exon[0], " ref_end: ", exon[1], " map_start: ", map_starts_l_l[j]," map_end: ", map_ends_l_l[j])
									if int(map_starts_l_l[j]) > int(exon[1]):
										break 
									read_region_l = [x for x in range(int(map_starts_l_l[j]),int(map_ends_l_l[j])+1)]
									read_region_set = set(read_region_l)

									exon_region_l = [x for x in range(int(exon[0]),int(exon[1])+1)]
									exon_region_set = set(exon_region_l)
	#								print("exon_region: ", exon_region_set)

									overlap = read_region_set & exon_region_set
									overlap_num = len(overlap)
									non_overlap = read_region_set ^ exon_region_set
									non_overlap_num = len(non_overlap)
									#non_overlap_in_map = read_region_set - exon_region_set
									#non_overlap_in_map_num = len(non_overlap_in_map)
									overlap_score = overlap_num - non_overlap_num
							
	##								print("overlap: ", overlap)
	##								print("overlap_num: ", overlap_num)
	##								print("non_overlap: ", non_overlap)
	##								print("non_overlap_num: ", non_overlap_num)
	#								print("non_overlap_in_map: ", non_overlap_in_map)
	#								print("non_overlap_in_map_num: ", non_overlap_in_map_num)
	##								print("overlap_score: ", overlap_score)
							

									if overlap_score > highest_overlap_score:#overlap領域(overlap数-非overlap数)がもっとも大きいもの
										highest_overlap_score = overlap_score
										exon_start_pos = int(exon[0])
										exon_end_pos = int(exon[1])
	##									print("exon_start_pos: ", exon_start_pos, sep="\t")
	##									print("highest_overlap_score: ", highest_overlap_score)
	#									print()
		
					
	##				print("highest_overlap_score", highest_overlap_score)
	##				print("exon_start_pos:",exon_start_pos)
	##				print("exon_end_pos:",exon_end_pos)
	##				print()	
					
					if not highest_overlap_score == -1000000:	
						#exon_len = (exon_end_pos - exon_start_pos) + 1
						"""
						print("exon_len: ", exon_len)
						print("read_len: ", read_len)
						print("cover_rate: ", cover_rate)
						print("")
						print("exon_start_pos: ", exon_start_pos, sep="\t")
						print("exon_end_pos: ", exon_end_pos, sep="\t")
						print("map_start_pos: ", map_starts_l_l[j], sep="\t")
						print("map_end_pos: ", map_ends_l_l[j], sep="\t")
						print(int(map_starts_l_l[j]) - (exon_start_pos+1), int(map_ends_l_l[j]) - exon_end_pos , sep="\t")
						print("")
						"""	
						
						if j_count == 1:
							if int(map_starts_l_l[j]) - (exon_start_pos) >= int(min_gap_1):#5'short
								if int(map_ends_l_l[j]) - exon_end_pos >= int(min_gap_2):#3'long
									read_condition += "slong,"
								else:
									read_condition += "short,"
							elif int(map_starts_l_l[j]) - (exon_start_pos) <= int("-" + min_gap_1):#5'long
								if int(map_ends_l_l[j]) - exon_end_pos <= int("-" + min_gap_2):#3'short
									read_condition += "slong,"
								else:
									read_condition += "long,"
							else:
								if int(map_ends_l_l[j]) - exon_end_pos >= int(min_gap_2):
									read_condition += "long,"
								elif int(map_ends_l_l[j]) - exon_end_pos <= int("-" + min_gap_2):
									read_condition += "short,"
								else:
									read_condition += "known,"
							
						elif j_count == len(map_starts_l_l):
							if int(map_ends_l_l[j]) - exon_end_pos >= int(min_gap_1):
								if int(map_starts_l_l[j]) - (exon_start_pos) >= int(min_gap_2):
									read_condition += "slong,"
								else:
									read_condition += "long,"
							elif int(map_ends_l_l[j]) - exon_end_pos <= int("-" + min_gap_1):
								if int(map_starts_l_l[j]) - (exon_start_pos) <= int("-" + min_gap_2):
									read_condition += "slong,"
								else:
									read_condition += "short,"
							else:
								if int(map_starts_l_l[j]) - (exon_start_pos) <= int("-" + min_gap_2):
									read_condition += "long,"
								elif int(map_starts_l_l[j]) - (exon_start_pos) >= int(min_gap_2):
									read_condition += "short,"
								else:
									read_condition += "known,"
						
						else:
							if int(map_starts_l_l[j]) - (exon_start_pos) <= int("-" + min_gap_2):
								if int(map_ends_l_l[j]) - exon_end_pos <= int("-" + min_gap_2):
									read_condition += "slong,"
								else:
									read_condition += "long,"
							elif int(map_starts_l_l[j]) - (exon_start_pos) >= int(min_gap_2):
								if int(map_ends_l_l[j]) - exon_end_pos >= int(min_gap_2):
									read_condition += "slong,"
								else:
									read_condition += "short,"
							else:
								if int(map_ends_l_l[j]) - exon_end_pos >= int(min_gap_2):
									read_condition += "long,"
								elif int(map_ends_l_l[j]) - exon_end_pos <= int("-" + min_gap_2):
									read_condition += "short,"
								else:
									read_condition += "known,"
						
					else:
	#					print("??")
						read_condition += "novel,"

	#				print(read_condition)
							
					if highest_overlap_score == -1000000:
						exon_position = "novel" + "\t" + "*" + "\t" + "*"
					else:
						exon_position = chr_l[i] + "\t" + str(exon_start_pos) + "\t" + str(exon_end_pos)
				
					annotated_exon_position.append(exon_position)
					highest_overlap_score = -1000000
	##			print("annotated_exon_position: ", annotated_exon_position)	
				"""
				#print(annotated_exon_position)
				>>> ['chrX\t119469660\t119470147', 'chrX\t119470372\t119470513', 'chrX\t119470900\t119471396']
				"""
	#			print(annotated_exon_position)
	#			annot = [exon_pos2RNAid_sort_dict[x] for x in annotated_exon_position]
				annot = [exon_pos2RNAid[x] for x in annotated_exon_position]
	#			for i in range(len(annot)):
	#				print(exon_pos2RNAid[i])
							
				annotated_exon_position = []
	#			print("annot: ", annot)
				"""
				print("annot: ", annot)
				>>> annot:  [[['COX6A1', '1']], [['AL021546.6', '2'], ['COX6A1', '2']], [['COX6A1', '5']]]
				"""
				
				annot_gene_num_dict = {}
				for i in range(len(annot)):
	#				print(annot[i])
					for j in range(len(annot[i])):
	#					print(annot[i][j])
						if annot[i][j][0] in annot_gene_num_dict:
							annot_gene_num_dict[annot[i][j][0]] += 1
						else:
							annot_gene_num_dict[annot[i][j][0]] = 1
				"""
				for k,v in sorted(annot_gene_num_dict.items(), key=lambda x: -x[1]):
					print(k,v,sep="\t")
				#COX6A1  3
				#AL021546.6      1
				"""
				annot_gene_num_dict_sort = sorted(annot_gene_num_dict.items(), key=lambda x: -x[1])
				
	#			print("annot_gene_num_dict_sort: ", annot_gene_num_dict_sort)
		
				gene_annot = OrderedDict()
				for i in range(len(annot)):
	#				print("annot[i]: ", annot[i], sep="\t")
					if annot[i][0] == "novel":
	#					print("novel")
						if not annot[i][0] in gene_annot:
							gene_annot[annot[i][0]] = annot[i][1]
						else:
							gene_annot[annot[i][0]] += "," + annot[i][1] 
					else:
						annot_dict = dict(annot[i])
	#					print("annot_dict: ", annot_dict)
						eva_annot = False
						for j in range(len(annot[i])):
	#						print("annot[i][j]: ", annot[i][j], sep="\t")
							if len(annot[i]) == 1:
								if not annot[i][j][0] in gene_annot:
									gene_annot[annot[i][j][0]] = annot[i][j][1]
								else:
									gene_annot[annot[i][j][0]] += ","+ annot[i][j][1]
							else:
								for k in range(len(annot_gene_num_dict_sort)):
									if eva_annot == False:
										if annot_gene_num_dict_sort[k][0] in annot_dict:
											eva_annot = True
											if not annot_gene_num_dict_sort[k][0] in gene_annot:
												gene_annot[annot_gene_num_dict_sort[k][0]] = annot_dict[annot_gene_num_dict_sort[k][0]]
											else:
												gene_annot[annot_gene_num_dict_sort[k][0]] += "," +  annot_dict[annot_gene_num_dict_sort[k][0]]
									

	#			print(gene_annot)
	#			print(len(gene_annot))
	#			if len(gene_annot) > 1:
	#				for i in range(len(gene_annot)):
						
	#			print(len(gene_annot))
				annotated_gene_len = len(gene_annot)

				read_id_l = line_l[0].split("/")
	#			print(read_id_l)

	#			print("len(map_starts_l)", len(map_starts_l), i_count)
				
				roop_count = 0
				for k,v in gene_annot.items():
					roop_count += 1
	#				print(k,v)
					v = v.split(",")
					v_str = ",".join(v)
	#				if annotated_gene_len == 1:
					if len(map_starts_l) == 1 and len(gene_annot) == 1: 
						print(line, k, v_str, read_condition, sep="\t")
						roop_count = 0
					else:
						annot_gene_name += k + "/" 
						annot_gene_exonNum += v_str + "/"
						read_conditions += read_condition + "/"
						read_condition = ""			
	#					if len(map_starts_l) == i_count: 
						if len(map_starts_l) == i_count and len(gene_annot) == roop_count: 
							print(line, annot_gene_name, annot_gene_exonNum, read_conditions, sep="\t")
							i_count = 0
							roop_count = 0
						
	#			print("")
	#	print("")

def annot_transcript(f, f_2):

	def del_terminal_index(l):
		if l[-1] == "":
			del l[-1]
		return l

	f1 = open(f)
	convert_gene_into_variant_dict = {}
	variant_strand_direction_dict = {}
	isoNum_dict = {}
	isoNum_uniq_dict = {}

	for line in f1:
		line = line.replace("\n","")
		line_l = line.split("\t")
		if not len(line_l) >= 2:
			continue 
		convert_gene_into_variant_dict[line_l[12],line_l[16]] = line_l[1]
		variant_strand_direction_dict[line_l[1]] = line_l[3]
		
		if not line_l[12] in isoNum_dict:
			exon_num_list = [line_l[16]]
			isoNum_dict[str(line_l[12])] = exon_num_list
		else:
			exon_num_list = [line_l[16]]
			isoNum_dict[str(line_l[12])] += exon_num_list

	for k ,v in isoNum_dict.items():
		uniq_v = list(set(v))
		isoNum_uniq_dict[k] = uniq_v

	isoNum_uniq_dict['novel'] = '*'
	isoNum_uniq_dict['FP1'] = '1'
	isoNum_uniq_dict['RP1'] = '1'
	isoNum_uniq_dict['FP2'] = '1'
	isoNum_uniq_dict['RP2'] = '1'
	isoNum_uniq_dict['P1'] = '1'
	isoNum_uniq_dict['P2'] = '1'
	#print(isoNum_uniq_dict['VEGFA'])
	#print(isoNum_dict['VEGFA'])

	f2 = open(f_2)

	for line in f2:
		line = line.replace("\n","")
		line_l = line.split("\t")
	
	#	print(line) 
	
		gene_l = line_l[11].split("/")
	#	print(gene_l)
		if len(gene_l) > 1:
			del gene_l[-1]
	
		exonNum_l = line_l[12].split("/")
	#	print(exonNum_l)
		if len(exonNum_l) > 1:
			del exonNum_l[-1]

		ref_exonNum_l = []
		for i in range(len(gene_l)):
			if gene_l[i] == "TRAV14":#refの遺伝子名に"/"が含まれる
				ref_exonNum_l.append(isoNum_uniq_dict[gene_l[i]+"/DV4"]) 
			elif gene_l[i] == "THRA1":
				ref_exonNum_l.append(isoNum_uniq_dict[gene_l[i]+"/BTR"])
			elif gene_l[i] in isoNum_uniq_dict:
				ref_exonNum_l.append(isoNum_uniq_dict[gene_l[i]])
		
	#	print("gene_l: ", gene_l)
	#	print("exonNum_l: ", exonNum_l)
	#	print("ref_exonNum_l: ", ref_exonNum_l)

	#	if line_l[9] in isoNum_uniq_dict[line_l[8]]:
		
		variant_name_s = ""
		variant_name_l = []
		TSSTES_s = ""
		TSSTES_l = []
		
	#	print("gene_l: ", gene_l)
		#gene_l:  ['FP1', 'MYL6B', 'MYL6']	

		for i in range(len(gene_l)):
			if gene_l[i] == "TRAV14":#refの遺伝子名に"/"が含まれる
				gene_l[i] = "TRAV14/DV4"
			elif gene_l[i] == "THRA1":
				gene_l[i] = "THRA1/BTR"
			elif not gene_l[i] in isoNum_uniq_dict:
				continue
			
			if gene_l[i] == 'FP1' or gene_l[i] == 'RP1' or gene_l[i] == 'FP2' or gene_l[i] == 'RP2' or gene_l[i] == 'P1' or gene_l[i] == 'P2':
				variant_name_l.append('PRIMER')
				TSSTES_l.append('*')
				continue	
			
			if gene_l[i] == 'novel':
				variant_name_l.append('UNANNOT')
				TSSTES_l.append('*')
				continue

			else:
				if exonNum_l[i] in ref_exonNum_l[i]:
					variant_name_l.append(convert_gene_into_variant_dict[gene_l[i], exonNum_l[i]])
					TSSTES_l.append('TSS+TES')
					continue
	#			print("ref_exonNum_l[i]: ", ref_exonNum_l[i])
				exonNum_l_l = exonNum_l[i].split(",")
				for j in range(len(ref_exonNum_l[i])):
	#				print("j: ", j)
	#				print("exonNum_l[i]: ", exonNum_l[i])
	#				print("gene_l[i]: ", gene_l[i])
	#				print("exonNum_l_l: ", exonNum_l_l)
	#				print("ref_exonNum_l[i]: ", ref_exonNum_l[i])
	#				print("ref_exonNum_l[i][j]: ", ref_exonNum_l[i][j])

					ref_exonNum_l_l = ref_exonNum_l[i][j].split(",")
	#				print("ref_exonNum_l_l: ", ref_exonNum_l_l)
					eva_match = False
					if set(exonNum_l_l) <= set(ref_exonNum_l_l):
						match_index_l = []
						for k in range(len(ref_exonNum_l_l)):
							if ref_exonNum_l_l[k] in exonNum_l_l:
								match_index_l.append(k)
	#					print("match_index_l: ", match_index_l, sep="\t")
						if len(match_index_l) == 1:		
							eva_match = True	
						else:
							for k in range(1, len(match_index_l)):
								if int(match_index_l[k]) == int(match_index_l[k-1]) + 1:
									eva_match = True
								else:
									eva_match = False
									break
	#				print(">>> eva_match: ", eva_match, sep="\t") 
		
	#				if set(exonNum_l_l) <= set(ref_exonNum_l_l):#181130変更
	#				if exonNum_s in ref_exonNum_s:#181130変更
					if eva_match == True:
	#					print(gene_l[i], ref_exonNum_l_l)
	#					print("exonNum_s: ", exonNum_s, sep="\t")	
	#					print("ref_econNum_s: ", ref_exonNum_s, sep="\t")
	#					print("exonNum_l_l: ", exonNum_l_l, sep="\t")
	#					print("ref_exonNum_l_l: ", ref_exonNum_l_l, sep="\t")
						if exonNum_l_l[0] == ref_exonNum_l_l[0] and exonNum_l_l[-1] != ref_exonNum_l_l[-1]:
							variant_name_s += convert_gene_into_variant_dict[gene_l[i], ref_exonNum_l[i][j]] + ","
							if variant_strand_direction_dict[convert_gene_into_variant_dict[gene_l[i], ref_exonNum_l[i][j]]] == "+":
								TSSTES_s += "TSS,"
							else:
								TSSTES_s += "TES,"
					
						elif exonNum_l_l[0] != ref_exonNum_l_l[0] and exonNum_l_l[-1] == ref_exonNum_l_l[-1]:
							variant_name_s += convert_gene_into_variant_dict[gene_l[i], ref_exonNum_l[i][j]] + ","
							if variant_strand_direction_dict[convert_gene_into_variant_dict[gene_l[i], ref_exonNum_l[i][j]]] == "+":
								TSSTES_s += "TES,"
							else:
								TSSTES_s += "TSS,"
						else:
							variant_name_s += convert_gene_into_variant_dict[gene_l[i], ref_exonNum_l[i][j]] + ","
							TSSTES_s += "None,"
					else:
	##					if j == 0:
						if j == len(ref_exonNum_l[i]) -1 and variant_name_s == "":
							variant_name_s += "Novel,"
							TSSTES_s += "*,"
	#						break
	#					else:
	#						break

				if not variant_name_s == "":
					variant_name_l.append(variant_name_s[:-1])			
					TSSTES_l.append(TSSTES_s[:-1])
					variant_name_s = ""
					TSSTES_s = ""
	#	print("variant_name_l: ", variant_name_l)
	#	print("TSSTES_l: ", TSSTES_l)
		
		print(line, "/".join(variant_name_l), "/".join(TSSTES_l), sep="\t")
	#	print("")

