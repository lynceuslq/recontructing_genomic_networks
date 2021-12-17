#!/bin/bash

fafile="GPD_sequences.fa"
torm_file_list="samprmlist.aa"
woringdir="/new_index"
tmptab=$woringdir/tmptabfile
BOWTIE2B="/bin/bowtie2-build"

echo -e "start to process full fa files"
#cat $fafile | tr "\n" "\t" | tr ">" "\n" | sed '/^[[:space:]]*$/d' > $tmptab

cat $torm_file_list | while read file
do

gmtoremove=$woringdir/$file.list
outfa=$woringdir/$file.fna
echo "start working on removing list $file at $(date)"

grep -v -w -f $gmtoremove $tmptab > $woringdir/${file}.tmptab
cat $woringdir/${file}.tmptab | tr " " "\t" | while read c1 c2 c3; do echo -e ">$c1\n$c3" ; done > $outfa
echo -e "ending working on $file at $(date)"

done
echo -e "genome extraction completed at $(date)"


cat $torm_file_list | while read file
do
outfa=$woringdir/$file.fna
outindex=$woringdir/$file
echo "start to generate index on $file at $(date)"
$BOWTIE2B $outfa $outindex --threads 8 --large-index

echo -e "cmpleted generating index on $file at $(date)"

done

echo -e "job complted at $(date)"
