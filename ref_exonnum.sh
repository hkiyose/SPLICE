#!/bin/bash
#module load python/3.6

REF=$1
REF_SEQ=$2
REF_DIR=$3

if [ ! -f $REF ]; then echo "$REF not found !!"; exit; fi
if [ ! -f $REF_SEQ ]; then echo "$REF_SEQ not found !!"; exit; fi
if [ ! -d $REF_DIR ]; then echo "$REF_DIR not found !!"; mkdir "$RED_DIR"; fi

echo "=========================================="
echo "Reference file; " $REF
echo "Reference genome sequence; " $REF_SEQ
echo "Directory of reference files; " $REF_DIR
echo "=========================================="

SRC=./src

python3 $SRC/add_exonnum.py $REF > $REF_DIR/`basename $REF`.exonnum

python3 $SRC/rm_conjoined.py $REF_DIR/`basename $REF`.exonnum > $REF_DIR/`basename $REF`.exonnum2

python3 $SRC/ref_convert.py $REF_DIR/`basename $REF`.exonnum2 > $REF_DIR/`basename $REF`.exonnum3

rm $REF_DIR/`basename $REF`.exonnum
rm $REF_DIR/`basename $REF`.exonnum2
mv $REF_DIR/`basename $REF`.exonnum3 $REF_DIR/`basename $REF`.exonnum

python3 $SRC/convert_fa.py $REF_SEQ > $REF_DIR/`basename $REF_SEQ`.liner

python3 $SRC/make_ref_seq_cp.py $REF_DIR/`basename $REF_SEQ`.liner $REF_DIR/`basename $REF`.exonnum > $REF_DIR/`basename $REF`.exonnum.fa

rm $REF_DIR/`basename $REF_SEQ`.liner

echo "FINISH!"
