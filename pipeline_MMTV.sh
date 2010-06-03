#!/bin/bash

# Blat to reference genome
# Build focused reference from blat
# Blat leftovers again to inserts

FNA=$1
BASE=$2
MKDB=$3

blat ~/mouseref/mouse_ref.fa $FNA -out=psl -minScore=10 $BASE.psl 2> $BASE.blatlog
exonerate -s 75 $FNA ~/marco/vecfiles/vec.fa > $BASE.exo &
exonerate -s 25 $FNA ~/marco/vecfiles/spMMTVbarcodes.fa > $BASE.bc.exo &
~/marco/pipeline/process_MMTV.pl $FNA $BASE $BASE $BASE $MKDB > $BASE.tsv 2> $BASE.log
~/marco/pipeline/makeRef.pl $BASE
blat $BASE.ref.fa $BASE.que.fa -minScore=10 -out=psl $BASE.2.psl
~/marco/pipeline/assign.pl $BASE.2.psl 200 $BASE >> $BASE.log
~/marco/pipeline/export_nogenes.pl $BASE > $BASE.export.txt
# ~/marco/pipeline/export.pl $BASE insert.db 1 > $BASE.export_squated.txt
