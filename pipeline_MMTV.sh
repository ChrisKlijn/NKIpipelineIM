#!/bin/bash

# Blat to reference genome
# Build focused reference from blat
# Blat leftovers again to inserts

FNA=$1
BASE=$2
MKDB=$3
BCFILE=$4

# Location of the pipeline and the reference files
# Change when applicable to a local install
reffile=/data/medoid-groups/im-mapping-pipeline/resource/mouse_ref.fa
pipelinedir=/data/medoid-groups/im-mapping-pipeline/NKIpipelineIM

blat $reffile $FNA -out=psl -minScore=10 $BASE.psl 2> $BASE.blatlog
exonerate -s 75 $FNA $pipelinedir/vecfiles/vec.fa > $BASE.exo &
exonerate -s 25 $FNA $BCFILE > $BASE.bc.exo &
$pipelinedir/process_MMTV.pl $FNA $BASE $BASE $BASE $MKDB > $BASE.tsv 2> $BASE.log
$pipelinedir/makeRef.pl $BASE
blat $BASE.ref.fa $BASE.que.fa -minScore=10 -out=psl $BASE.2.psl
$pipelinedir/assign.pl $BASE.2.psl 200 $BASE >> $BASE.log
$pipelinedir/export_nogenes.pl $BASE > $BASE.export.txt
# $pipelinedir/export.pl $BASE insert.db 1 > $BASE.export_squated.txt
