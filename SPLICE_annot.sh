#!/bin/bash
module load python/3.6

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

echo ">>> Base quality filtering (1/20)"
echo "python $SRC/bq_filter.py $FASTQ $BQ_FILT > $OUT/`basename $FASTQ`.bqf"
python $SRC/bq_filter.py $FASTQ $BQ_FILT > $OUT/`basename $FASTQ`.bqf
echo ""

echo ">>> Mapping to reference genome (2/20)"
echo "$MINIMAP2 -ax splice --cs $REF_GENOME_FA $OUT/`basename $FASTQ`.bqf > $OUT/`basename $FASTQ`.bqf.sam"
$MINIMAP2 -ax splice --cs $REF_GENOME_FA $OUT/`basename $FASTQ`.bqf > $OUT/`basename $FASTQ`.bqf.sam
echo ""

echo ">>> Convert the SAM(genome) format (3/20)"
echo "python $SRC/sam_convert.py $OUT/`basename $FASTQ`.bqf.sam | python $SRC/sam_convert2.py /dev/stdin $OUT/`basename $FASTQ`.bqf 0.3 > $OUT/`basename $FASTQ`.bqf.sam.conv"
python $SRC/sam_convert.py $OUT/`basename $FASTQ`.bqf.sam | python $SRC/sam_convert2.py /dev/stdin $OUT/`basename $FASTQ`.bqf 0.3 > $OUT/`basename $FASTQ`.bqf.sam.conv
echo ""

echo ">>> Re-mapping the softclip region (4/20)"
echo "python $SRC/softclip_fa.py $OUT/`basename $FASTQ`.bqf.sam.conv $OUT/`basename $FASTQ`.bqf $MIN_SC_LEN > $OUT/`basename $FASTQ`.bqf.sam.conv.fa"
python $SRC/softclip_fa.py $OUT/`basename $FASTQ`.bqf.sam.conv $OUT/`basename $FASTQ`.bqf $MIN_SC_LEN > $OUT/`basename $FASTQ`.bqf.sam.conv.fa

echo "$MINIMAP2 -ax splice --cs $REF_GENOME_FA $OUT/`basename $FASTQ`.bqf.sam.conv.fa > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam"
$MINIMAP2 -ax splice --cs $REF_GENOME_FA $OUT/`basename $FASTQ`.bqf.sam.conv.fa > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam

rm $OUT/`basename $FASTQ`.bqf.sam.conv.fa
echo ""

echo ">>> Integrating mapped softclip regions into converted SAM (5/20)"
echo "python $SRC/sam_convert.py $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv"
python $SRC/sam_convert.py $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv

echo "cat $OUT/`basename $FASTQ`.bqf.sam.conv >> $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv"
cat $OUT/`basename $FASTQ`.bqf.sam.conv >> $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv

echo "sort -k1 $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv | python $SRC/sam_convert3.py /dev/stdin $OUT/`basename $FASTQ`.bqf 0.3 > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2"
sort -k1 $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv | python $SRC/sam_convert3.py /dev/stdin $OUT/`basename $FASTQ`.bqf 0.3 > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2

rm $OUT/`basename $FASTQ`.bqf.sam.conv
rm $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv
echo ""

echo ">>> Gene annotation (6/20)"
echo "python $SRC/annot_exonnum.py $REF_TRANSCRIPT $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2 $MQ_FILT 1 1 > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot"
python $SRC/annot_exonnum.py $REF_TRANSCRIPT $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2 $MQ_FILT 1 1 > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot
echo ""

echo ">>> Splicing variant annotation (7/20)"
echo "python $SRC/annot_transcript.py $REF_TRANSCRIPT $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2"
python $SRC/annot_transcript.py $REF_TRANSCRIPT $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2
echo ""

echo ">>> Mapping to reference transcriptome (8/20)"
echo "$MINIMAP2 -a --cs $REF_TRANSCRIPT_FA $OUT/`basename $FASTQ`.bqf > $OUT/`basename $FASTQ`.bqf.cDNA.sam"
$MINIMAP2 -a --cs $REF_TRANSCRIPT_FA $OUT/`basename $FASTQ`.bqf > $OUT/`basename $FASTQ`.bqf.cDNA.sam
echo ""

echo ">>> Add information on the number of alignment matches (Converted SAM(Genome)/SAM(Transcriptome)) (9/20)"
echo "python $SRC/match_rate.py $OUT/`basename $FASTQ`.bqf.cDNA.sam > $OUT/`basename $FASTQ`.bqf.cDNA.sam.match_num"
python $SRC/match_rate.py $OUT/`basename $FASTQ`.bqf.cDNA.sam > $OUT/`basename $FASTQ`.bqf.cDNA.sam.match_num

echo "python $SRC/calc_error.py $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2 > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj"
python $SRC/calc_error.py $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2 > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj

echo "python $SRC/match_rate2.py $OUT/`basename $FASTQ`.bqf.cDNA.sam.match_num $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2"
python $SRC/match_rate2.py $OUT/`basename $FASTQ`.bqf.cDNA.sam.match_num $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2

rm $OUT/`basename $FASTQ`.bqf.cDNA.sam.match_num
echo ""

echo ">>> Table (10/20)"
echo "python $SRC/format_change.py $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2 > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table"
python $SRC/format_change.py $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2 > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table
echo ""

echo ">>> Merge at position inside the TSS side exon and inside the TES side exon (11/20)"
echo "sort -k3 $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table | python $SRC/merge_transcript.py $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2 /dev/stdin > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv"
sort -k3 $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table | python $SRC/merge_transcript.py $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2 /dev/stdin > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv

echo "python $SRC/merge_transcript2.py $REF_TRANSCRIPT $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2 > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap"
python $SRC/merge_transcript2.py $REF_TRANSCRIPT $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2 > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap

echo "python $SRC/merge_transcript3.py $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg"
python $SRC/merge_transcript3.py $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg
echo ""

echo ">>> Detect novel exon combination candidates (12/20)"
echo "sort -nr $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg | python $SRC/det_novel.py $REF_TRANSCRIPT /dev/stdin mid_comb > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg.mid_comb"
sort -nr $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg | python $SRC/det_novel.py $REF_TRANSCRIPT /dev/stdin mid_comb > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg.mid_comb
echo ""

echo ">>> Detect novel exon length candidates (13/20)"
echo "sort -nr $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg | python $SRC/det_novel.py $REF_TRANSCRIPT /dev/stdin mid_sl > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg.mid_sl"
sort -nr $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg | python $SRC/det_novel.py $REF_TRANSCRIPT /dev/stdin mid_sl > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg.mid_sl
echo ""

echo ">>> Detect novel TES candidates (14/20)"
echo "python $SRC/det_novel_tes.py $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2"
python $SRC/det_novel_tes.py $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2

echo "sort -nr $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2 | python $SRC/det_novel_tes2.py $REF_TRANSCRIPT /dev/stdin > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known"
sort -nr $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2 | python $SRC/det_novel_tes2.py $REF_TRANSCRIPT /dev/stdin > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known

echo "sort -k3 $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known | python $SRC/det_novel_tes3.py $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2 /dev/stdin > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known2"
sort -k3 $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known | python $SRC/det_novel_tes3.py $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2 /dev/stdin > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known2

echo "python $SRC/det_novel_tes4.py $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known2 > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known3"
python $SRC/det_novel_tes4.py $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known2 > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known3

echo "python $SRC/det_novel_tes5.py $REF_TRANSCRIPT $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known3 > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known4"
python $SRC/det_novel_tes5.py $REF_TRANSCRIPT $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known3 > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known4

echo "python $SRC/det_novel_tes6.py $REF_TRANSCRIPT $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2 $OUT/`basename $FASTQ`.bqf 15 0.8 $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known4 > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known4.pA"
python $SRC/det_novel_tes6.py $REF_TRANSCRIPT $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2 $OUT/`basename $FASTQ`.bqf 15 0.8 $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known4 > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known4.pA

echo "sort -nr -k2 $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known4.pA | python $SRC/det_novel_tes7.py $REF_TRANSCRIPT /dev/stdin $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2 > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known4.pA.TES"
sort -nr -k2 $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known4.pA | python $SRC/det_novel_tes7.py $REF_TRANSCRIPT /dev/stdin $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2 > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known4.pA.TES
echo ""

echo ">>> Convert the SAM(transcriptome) format (15/20)"
echo "python $SRC/sam_convert.py $OUT/`basename $FASTQ`.bqf.cDNA.sam | python $SRC/sam_convert2.py /dev/stdin $OUT/`basename $FASTQ`.bqf 0.3 > $OUT/`basename $FASTQ`.bqf.cDNA.sam.conv"
python $SRC/sam_convert.py $OUT/`basename $FASTQ`.bqf.cDNA.sam | python $SRC/sam_convert2.py /dev/stdin $OUT/`basename $FASTQ`.bqf 0.3 > $OUT/`basename $FASTQ`.bqf.cDNA.sam.conv
echo ""

echo ">>> Detect fusion gene candidates (16/20)"
echo "python $SRC/det_fusion.py $REF_TRANSCRIPT $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2 $OUT/`basename $FASTQ`.bqf.cDNA.sam.conv $OUT/`basename $FASTQ`.bqf 15 0.8 > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.fusion"
python $SRC/det_fusion.py $REF_TRANSCRIPT $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2 $OUT/`basename $FASTQ`.bqf.cDNA.sam.conv $OUT/`basename $FASTQ`.bqf 15 0.8 > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.fusion

echo "python $SRC/det_fusion2.py $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.fusion $MIN_FUSION_DIST | python $SRC/det_fusion3.py /dev/stdin $MAX_FUSION_BP_MERGE $MIN_FUSION_READ $MIN_FUSION_FREQ | sort -nr > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.fusion2"
python $SRC/det_fusion2.py $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.fusion $MIN_FUSION_DIST | python $SRC/det_fusion3.py /dev/stdin $MAX_FUSION_BP_MERGE $MIN_FUSION_READ $MIN_FUSION_FREQ | sort -nr > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.fusion2
echo ""

echo ">>> Add information on the number of alignment matches (Converted SAM(Transcriptome)) (17/20)"
echo "python $SRC/sam_convert4.py $OUT/`basename $FASTQ`.bqf.cDNA.sam.conv > $OUT/`basename $FASTQ`.bqf.cDNA.sam.conv.map_rate"
python $SRC/sam_convert4.py $OUT/`basename $FASTQ`.bqf.cDNA.sam.conv > $OUT/`basename $FASTQ`.bqf.cDNA.sam.conv.map_rate

rm $OUT/`basename $FASTQ`.bqf.cDNA.sam.conv
echo ""

echo ">>> Detect novel position candidates (18/20)"
echo "python $SRC/det_novel2.py $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2 $MQ_FILT_NOVEL_EXON $OUT/`basename $FASTQ`.bqf.cDNA.sam.conv.map_rate 0.8 60 $OUT/`basename $FASTQ`.bqf 15 0.8 > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.novel_exon"
python $SRC/det_novel2.py $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2 $MQ_FILT_NOVEL_EXON $OUT/`basename $FASTQ`.bqf.cDNA.sam.conv.map_rate 0.8 60 $OUT/`basename $FASTQ`.bqf 15 0.8 > $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.novel_exon

rm $OUT/`basename $FASTQ`.bqf.cDNA.sam.conv.map_rate
echo ""

echo ">>> Expression table (19/20)"
echo "python $SRC/get_table.py \
$OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg.mid_comb \
$OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg.mid_sl \
$OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known4.pA.TES \
$OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known4 \
$OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2 \
> $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.expression"

python $SRC/get_table.py \
$OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg.mid_comb \
$OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg.mid_sl \
$OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known4.pA.TES \
$OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.conv.gap.mrg2.known4 \
$OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2 \
> $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.expression
echo ""

echo ">>> Analysis report (20/20)"
echo "python $SRC/get_report.py $FASTQ $OUT/`basename $FASTQ`.bqf $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2 > $OUT/`basename $FASTQ`.bqf.stat"
python $SRC/get_report.py $FASTQ $OUT/`basename $FASTQ`.bqf $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2 > $OUT/`basename $FASTQ`.stat
echo ""

echo "mv $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.fusion2 $OUT/`basename $FASTQ`.fusion"
mv $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.fusion2 $OUT/`basename $FASTQ`.fusion
echo ""

echo "mv $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.novel_exon $OUT/`basename $FASTQ`.novel_exon"
mv $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.novel_exon $OUT/`basename $FASTQ`.novel_exon
echo ""

echo "mv $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.expression $OUT/`basename $FASTQ`.annot"
mv $OUT/`basename $FASTQ`.bqf.sam.conv.fa.sam.conv2.annot2.sj2.table.expression $OUT/`basename $FASTQ`.annot
echo ""

rm -r $OUT/`basename $FASTQ`.bqf*


echo "FINISH"
