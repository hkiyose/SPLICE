#!/bin/bash

IN=$1
OUT=$2

if [ ! -d $IN ]; then echo "$IN not found !!"; exit; fi
if [ ! -d $OUT ]; then echo "$OUT not found !!"; mkdir $OUT; fi

echo "=========================================="
echo "IN; " $IN
echo "OUT; " $OUT
echo "=========================================="

SRC=./src
CONFIG=./configure
. $CONFIG

#echo ">>> Merging tables (1/14)"
python3 $SRC/merge_table.py $IN > $OUT/expression.table

#echo ">>> Filtering out low-frequency Novel splicing variants (2/14)"
python3 $SRC/freq_filter.py $OUT/expression.table $MIN_READ_NUM $MIN_READ_FREQ > $OUT/expression.table2

#echo ">>> Change the format, evaluate the splicing junction error rate and remove unreliable novel splicing variants (3/14)"
python3 $SRC/format_change2.py $SRC/gene_type/genecode.txt $SRC/gene_type/refseq.txt $REF_TRANSCRIPT $OUT/expression.table2 30 $ERR_RATE_FILT $MIN_NOVEL_LEN_GAP $RANGE_SJ_EVA $RANGE_SJ_EVA > $OUT/expression.table3

#echo ">>> Merge each candidate by category (4/14)"
#Known candidates
python3 $SRC/merge_known.py $OUT/expression.table3 > $OUT/expression.table3.known

#Novel exon length candidates
python3 $SRC/merge_len_comb.py $REF_TRANSCRIPT $OUT/expression.table3 LENGTH | sort -k 4,4 -k 7nr,7 > $OUT/expression.table3.length

python3 $SRC/merge_len_comb2.py $OUT/expression.table3.length > $OUT/expression.table3.length2

python3 $SRC/merge_len_comb3.py $OUT/expression.table3.length2 > $OUT/expression.table3.length3

#Novel exon combination candidates
python3 $SRC/merge_len_comb.py $REF_TRANSCRIPT $OUT/expression.table3 COMBINATION | sort -k 4,4 -k 7nr,7 > $OUT/expression.table3.combination

python3 $SRC/merge_len_comb2.py $OUT/expression.table3.combination > $OUT/expression.table3.combination2

python3 $SRC/merge_len_comb3.py $OUT/expression.table3.combination2 > $OUT/expression.table3.combination3

#echo ">>> Merging novel unannotated exon tables (5/14)"
python3 $SRC/merge_novel_exon.py $IN > $OUT/novel_exon.table

python3 $SRC/merge_novel_exon2.py $OUT/novel_exon.table > $OUT/novel_exon.table2

python3 $SRC/merge_novel_exon3.py $OUT/novel_exon.table2 > $OUT/novel_exon.table3

python3 $SRC/merge_novel_exon4.py $OUT/novel_exon.table3 $OUT/novel_exon.table > $OUT/novel_exon.table4

python3 $SRC/merge_novel_exon5.py $OUT/novel_exon.table4 > $OUT/novel_exon.table4.convert

#echo ">>> Filtering out low-frequency Novel unannot exon candidates (6/14)"
python3 $SRC/freq_filter2.py $OUT/expression.table $OUT/novel_exon.table4.convert > $OUT/novel_exon.table5

#echo ">>> Extracting Novel unannot exon candidates detected together with the gene region (7/14)"
python3 $SRC/novel_exon_filter.py $SRC/gene_type/genecode.txt $SRC/gene_type/refseq.txt $OUT/novel_exon.table5 > $OUT/novel_exon.table6

#echo ">>> Removal of Fusion gene candidates containing Novel exon (8/14)"
sort -k4 $OUT/novel_exon.table6 | python3 $SRC/novel_exon_filter2.py $REF_TRANSCRIPT /dev/stdin | sort -k 4,4 -k 7nr,7 > $OUT/novel_exon.table7

#echo ">>> Merging novel unannot exon candidates (9/14)"
python3 $SRC/merge_novel_exon6.py $OUT/novel_exon.table7 > $OUT/novel_exon.table8

python3 $SRC/merge_len_comb3.py $OUT/novel_exon.table8 > $OUT/novel_exon.table9

#echo ">>> Conversion of novel exon location information to the most frequent location, filtering out unreliable candidates (10/14)"
python3 $SRC/novel_exon_sort.py $OUT/novel_exon.table > $OUT/novel_exon.table.sort

python3 $SRC/novel_exon_convert.py $OUT/novel_exon.table9 > $OUT/novel_exon.table10

python3 $SRC/novel_exon_convert2.py $OUT/novel_exon.table9 $OUT/novel_exon.table.sort $REF_TRANSCRIPT $OUT/novel_exon.table10 > $OUT/novel_exon.table11

python3 $SRC/novel_exon_convert3.py $REF_TRANSCRIPT $OUT/novel_exon.table11 $MIN_NOEL_EXON_LEN $MIN_NOEL_EXON_LEN $OUT/novel_exon.table9 > $OUT/novel_exon.table12

#echo ">>> Integrate each table (11/14)"
cp $OUT/expression.table3.known $OUT/expression.table3.merge

cat $OUT/novel_exon.table12 >> $OUT/expression.table3.merge

cat $OUT/expression.table3.combination3 >> $OUT/expression.table3.merge

cat $OUT/expression.table3.length3 >> $OUT/expression.table3.merge 

#echo ">>> Add information such as splicing variant frequency and total read count. Exclude low-frequency candidates (12/14)"
python3 $SRC/format_change3.py $OUT/expression.table $OUT/expression.table3.merge > $OUT/expression.table3.merge2

python3 $SRC/format_change4.py $OUT/expression.table3.merge2 $MIN_READ_FREQ $MIN_READ_NUM > $OUT/expression.table3.merge3

#echo ">>> Integrate 3'UTR information (13/14)"
python3 $SRC/merge_3utr.py $OUT/expression.table3.merge3 > $OUT/expression.table3.merge4

#echo ">>> Modify the coding type and generate the final file (14/14)"
python3 $SRC/format_change5.py $SRC/gene_type/genecode.txt $SRC/gene_type/refseq.txt $REF_TRANSCRIPT $OUT/expression.table3.merge4 > $OUT/expression.table3.merge5

python3 $SRC/format_change6.py $OUT/expression.table3.merge5 > $OUT/expression.tsv

rm $OUT/expression.table*
rm $OUT/novel_exon.table*

echo "OUTPUT: " $OUT/expression.tsv 
