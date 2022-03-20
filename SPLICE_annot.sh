#!/bin/bash
#module load python/3.6

FASTQ=$1
OUT=$2

if [ ! -f $FASTQ ]; then echo "$FASTQ not found !!"; exit; fi
if [ ! -d $OUT ]; then echo "$OUT not found !!"; mkdir $OUT; fi

echo "=========================================="
echo "FASTQ; " $FASTQ
echo "OUT; " $OUT
echo "=========================================="

SRC=./src 
CONFIG=./configure
. $CONFIG

#echo ">>> Base quality filtering (1/20)"
python3 $SRC/bq_filter.py $FASTQ $BQ_FILT > $OUT/`basename $FASTQ`.bqf

#echo ">>> Mapping to reference genome (2/20)"
$MINIMAP2 -ax splice --cs $REF_GENOME_FA $OUT/`basename $FASTQ`.bqf > $OUT/`basename $FASTQ`.bqf.sam

#echo ">>> Convert the SAM(genome) format (3/20)"
python3 $SRC/sam_convert.py $OUT/`basename $FASTQ`.bqf.sam | python3 $SRC/sam_convert2.py /dev/stdin $OUT/`basename $FASTQ`.bqf 0.3 > $OUT/`basename $FASTQ`.bqf.sam.conv

#echo ">>> Re-mapping the softclip region (4/20)"
python3 $SRC/softclip_fa.py $OUT/`basename $FASTQ`.bqf.sam.conv $OUT/`basename $FASTQ`.bqf $MIN_SC_LEN > $OUT/`basename $FASTQ`.bqf.sam.conv.fa

$MINIMAP2 -ax splice --cs $REF_GENOME_FA $OUT/`basename $FASTQ`.bqf.sam.conv.fa > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam

rm $OUT/`basename $FASTQ`.bqf.sam.conv.fa

#echo ">>> Integrating mapped softclip regions into converted SAM (5/20)"
python3 $SRC/sam_convert.py $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv

cat $OUT/`basename $FASTQ`.bqf.sam.conv >> $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv

sort -k1 $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv | python3 $SRC/sam_convert3.py /dev/stdin $OUT/`basename $FASTQ`.bqf 0.3 > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2

rm $OUT/`basename $FASTQ`.bqf.sam.conv
rm $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv

#echo ">>> Gene annotation (6/20)"
python3 $SRC/annot_exonnum.py $REF_TRANSCRIPT $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2 $MQ_FILT 1 1 > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot

#echo ">>> Splicing variant annotation (7/20)"
python3 $SRC/annot_transcript.py $REF_TRANSCRIPT $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2

#echo ">>> Mapping to reference transcriptome (8/20)"
$MINIMAP2 -a --cs $REF_TRANSCRIPT_FA $OUT/`basename $FASTQ`.bqf > $OUT/`basename $FASTQ`.bqf.cDNA.sam

#echo ">>> Add information on the number of alignment matches (Converted SAM(Genome)/SAM(Transcriptome)) (9/20)"
python3 $SRC/match_rate.py $OUT/`basename $FASTQ`.bqf.cDNA.sam > $OUT/`basename $FASTQ`.bqf.cDNA.sam.match_num

python3 $SRC/calc_error.py $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2 > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj

python3 $SRC/match_rate2.py $OUT/`basename $FASTQ`.bqf.cDNA.sam.match_num $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2

rm $OUT/`basename $FASTQ`.bqf.cDNA.sam.match_num

#echo ">>> Table (10/20)"
python3 $SRC/format_change.py $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2 > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table

#echo ">>> Merge at position inside the TSS side exon and inside the TES side exon (11/20)"
sort -k3 $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table | python3 $SRC/merge_transcript.py $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2 /dev/stdin > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv

python3 $SRC/merge_transcript2.py $REF_TRANSCRIPT $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2 > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap

python3 $SRC/merge_transcript3.py $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg

#echo ">>> Detect novel exon combination candidates (12/20)"
sort -nr $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg | python3 $SRC/det_novel.py $REF_TRANSCRIPT /dev/stdin mid_comb > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg.mid_comb

#echo ">>> Detect novel exon length candidates (13/20)"
sort -nr $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg | python3 $SRC/det_novel.py $REF_TRANSCRIPT /dev/stdin mid_sl > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg.mid_sl

#echo ">>> Detect novel TES candidates (14/20)"
python3 $SRC/det_novel_tes.py $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2

sort -nr $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2 | python3 $SRC/det_novel_tes2.py $REF_TRANSCRIPT /dev/stdin > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known

sort -k3 $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known | python3 $SRC/det_novel_tes3.py $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2 /dev/stdin > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known2

python3 $SRC/det_novel_tes4.py $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known2 > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known3

python3 $SRC/det_novel_tes5.py $REF_TRANSCRIPT $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known3 > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known4

python3 $SRC/det_novel_tes6.py $REF_TRANSCRIPT $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2 $OUT/`basename $FASTQ`.bqf 15 0.8 $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known4 > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known4.pA

sort -nr -k2 $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known4.pA | python3 $SRC/det_novel_tes7.py $REF_TRANSCRIPT /dev/stdin $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2 > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known4.pA.TES

#echo ">>> Convert the SAM(transcriptome) format (15/20)"
python3 $SRC/sam_convert.py $OUT/`basename $FASTQ`.bqf.cDNA.sam | python3 $SRC/sam_convert2.py /dev/stdin $OUT/`basename $FASTQ`.bqf 0.3 > $OUT/`basename $FASTQ`.bqf.cDNA.sam.conv

#echo ">>> Detect fusion gene candidates (16/20)"
python3 $SRC/det_fusion.py $REF_TRANSCRIPT $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2 $OUT/`basename $FASTQ`.bqf.cDNA.sam.conv $OUT/`basename $FASTQ`.bqf 15 0.8 > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.fusion

python3 $SRC/det_fusion2.py $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.fusion $MIN_FUSION_DIST | python3 $SRC/det_fusion3.py /dev/stdin $MAX_FUSION_BP_MERGE $MIN_FUSION_READ $MIN_FUSION_FREQ | sort -nr > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.fusion2

#echo ">>> Add information on the number of alignment matches (Converted SAM(Transcriptome)) (17/20)"
python3 $SRC/sam_convert4.py $OUT/`basename $FASTQ`.bqf.cDNA.sam.conv > $OUT/`basename $FASTQ`.bqf.cDNA.sam.conv.map_rate

rm $OUT/`basename $FASTQ`.bqf.cDNA.sam.conv

#echo ">>> Detect novel position candidates (18/20)"
python3 $SRC/det_novel2.py $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2 $MQ_FILT_NOVEL_EXON $OUT/`basename $FASTQ`.bqf.cDNA.sam.conv.map_rate 0.8 60 $OUT/`basename $FASTQ`.bqf 15 0.8 > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.novel_exon

rm $OUT/`basename $FASTQ`.bqf.cDNA.sam.conv.map_rate

#echo ">>> Expression table (19/20)"
python3 $SRC/get_table.py \
$OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg.mid_comb \
$OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg.mid_sl \
$OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known4.pA.TES \
$OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known4 \
$OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2 \
> $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.expression

#echo ">>> Analysis report (20/20)"
python3 $SRC/get_report.py $FASTQ $OUT/`basename $FASTQ`.bqf $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2 > $OUT/`basename $FASTQ`.stat

mv $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.fusion2 $OUT/`basename $FASTQ`.fusion

mv $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.novel_exon $OUT/`basename $FASTQ`.novel_exon

mv $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.expression $OUT/`basename $FASTQ`.annot

rm -r $OUT/`basename $FASTQ`.bqf*

echo "FINISH"
