#!/bin/bash

# Blat to reference genome
# Build focused reference from blat
# Blat leftovers again to inserts

FNA=$1
BASE=$2
MKDB=$3

# blat ~/mouseref/mouse_ref.fa $FNA -out=psl -minScore=10 $BASE.psl 2> $BASE.blatlog
# exonerate -s 75 $FNA ~/marco/vecfiles/vec.fa > $BASE.exo &
exonerate -s 49 --exhaustive yes $FNA ~/marco/vecfiles/MULVbc.fa > $BASE.bc.exo &
~/shearing/pipeline/process_MULV.pl $FNA $BASE $BASE $BASE $MKDB > $BASE.tsv 2> $BASE.log
~/shearing/pipeline/makeRef.pl $BASE
blat $BASE.ref.fa $BASE.que.fa -minScore=10 -out=psl $BASE.2.psl
~/shearing/pipeline/assign.pl $BASE.2.psl 200 $BASE >> $BASE.log
~/shearing/pipeline/export_nogenes.pl $BASE > $BASE.export.txt
#~/marco/pipeline/export.pl $BASE insert.db 1 > $BASE.export_squated.txt
