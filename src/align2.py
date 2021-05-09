
def SW(seq1, seq2, match, mismatch, gap):

    score = []
    pointer = []
    for i in range(len(seq1)):
        tmp = []
        tmp2 = []
        for j in range(len(seq2)):
            tmp.append(0)
            tmp2.append("")

        score.append(tmp)
        pointer.append(tmp2)

    max_score = 0
    max_i = 0
    max_j = 0
    for i in range(len(seq1)):
        for j in range(len(seq2)):
            if i - 1 >= 0 and j - 1 >= 0:
                (diagonal_score, left_score, up_score) = (score[i - 1][j - 1], score[i][j - 1] + gap, score[i - 1][j] + gap)
            else:
                (diagonal_score, left_score, up_score) = (0, 0 + gap, 0 + gap)

            if seq1[i] == seq2[j]:
                diagonal_score += match
            else:
                diagonal_score += mismatch

            if diagonal_score > left_score and diagonal_score > up_score:
                pointer[i][j] = "d"
                if diagonal_score > 0:
                    score[i][j] = diagonal_score
                else:
                    score[i][j] = 0
            elif left_score > diagonal_score and left_score > up_score:
                pointer[i][j] = "l"
                if left_score > 0:
                    score[i][j] = left_score
                else:
                    score[i][j] = 0
            else:
                pointer[i][j] = "u"
                score[i][j] = up_score
                if up_score > 0:
                    score[i][j] = up_score
                else:
                    score[i][j] = 0

            if max_score < score[i][j]:
                max_score = score[i][j]
                max_i = i
                max_j = j

#    seq_list = [" "]
#    for tmp in list(seq2):
#        seq_list.append(tmp+" ")
#    print(seq_list)
#    for i in range(len(seq1)):
#        out = [seq1[i]]
#        for j in range(len(seq2)):
#             out.append(str(score[i][j]) + pointer[i][j])
#        print(out)

    seq1_alignd = ""
    seq2_alignd = ""
    i = max_i
    j = max_j
    i2 = 0
    j2 = 0	
    while i >=0 and j >= 0 and score[i][j] > 0:
        if pointer[i][j] == "d":
            seq1_alignd += seq1[i]
            seq2_alignd += seq2[j]
            i2 += 1
            j2 += 1
            i = i - 1
            j = j - 1
        elif pointer[i][j] == "l":
            seq1_alignd += "-"
            seq2_alignd += seq2[j]
            j2 += 1
            i = i
            j = j - 1
        elif pointer[i][j] == "u":
            seq1_alignd += seq1[i]
            seq2_alignd += "-"
            i2 += 1
            i = i - 1
            j = j

    return seq1_alignd[::-1], seq2_alignd[::-1], max_i - i2 + 1, max_i, max_j - j2 + 1, max_j

if __name__ == '__main__':
	seq1 = "AACTT"
	seq2 = "AATT"
	match = 1
	mismatch = -1
	gap = -1

	res = []
	res = SW(seq1, seq2, match, mismatch, gap)

	print(seq1)
	print(seq2)

	print(res)
	print("SW alignment")
	print(res[0])
	print(res[1])

