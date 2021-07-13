#!/bin/bash

export PATH=/PATH/TO/BLAST/ncbi-blast-2.10.1+/bin:$PATH

MCXLOAD="/PATH/TO/mcxload"
MCL="/PATH/TO/mcl"
MCXDUMP="/PATH/TO/mcxdump"


workingdir=""
INPUTdir="/PATH/TO/gut_phage_database"
phagelist="/PATH/TO/splitbyhosttaxa.withhost.txt"#phage genome accessions and host infomation, separated by a tab

echo -e "start to prepare fasta files at $(date)"

cut -f1 $phagelist | sort | uniq | while read acc; do grep -w "$acc" -A1  $INPUTdir/GPD_sequences.fa >> $workingdir/tmpin.fa  ; done

echo -e "finishing merging $(grep -c ">" $workingdir/tmpin.fa) sequences from hosts at $(date)"

INPUTfna=$workingdir/tmpin.fa
OUTPUTtab=$workingdir/tmpin.outfmt6
OUTPUTabc=$workingdir/tmpin.abc

echo -e "start to generate database at $(date)"

makeblastdb -in $INPUTfna -out $workingdir/all2all_tmp -dbtype nucl

echo -e "start to blast at $(date)"

blastn -query $INPUTfna -db $workingdir/all2all_tmp -evalue 10 -perc_identity 99 -num_threads 8 -outfmt "6 qacc sacc qlen slen length pident evalue bitscore mismatch qstart qend sstart send qseq sseq" -out $workingdir/tmp.out

echo -e "start to filter results at $(date)"

awk '$1 != $2' $workingdir/tmp.out | awk '$5 >= 500' > $OUTPUTtab

cut -f 1,2,7 $OUTPUTtab | sort | uniq > $OUTPUTabc

cat $OUTPUTtab | cut -f1,2,5 | sort -u > $workingdir/blastn.length.abc

echo -e "finishing slicing blastoutput at $(date)"

cut -f1,2 $workingdir/blastn.length.abc | sort -u | while read c1 c2 ; do awk -v a="$c1" -v b="$c2" '$1 == a && $2 ==b' $workingdir/blastn.length.abc | sort -k3,3 -n -r | head -1 ; done > $workingdir/blastn.length.sorted.top.abc

echo -e "finishing sorting out top alignment length at $(date)"



echo -e "start to load abc file at $(date)"

$MCXLOAD  -abc $workingdir/tmpin.abc --stream-mirror  --stream-neg-log10 -stream-tf 'ceil(200)' -write-tab $workingdir/tmpin.tab -o $workingdir/tmpin.mci

echo -e "start to generate clusters at $(date)"

$MCL $workingdir/tmpin.mci -I 1.4

echo -e "start to make tabular results at $(date)"

$MCXDUMP -icl $workingdir/out.tmpin.mci.I14 -tabr  $workingdir/tmpin.tab -o  $workingdir/clusters.mci.I14

echo -e "MCL completed at $(date)"

echo -e "job completed at $(date)"
