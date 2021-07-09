#!/bin/bash

export PATH=/hwfssz5/ST_INFECTION/GlobalDatabase/user/liqian6/tools/ncbi-blast-2.10.1+/bin:$PATH

MCXLOAD="/hwfssz5/ST_INFECTION/GlobalDatabase/share/software/Miniconda3/envs.multi-user/vcontact2/bin/mcxload"
MCL="/hwfssz5/ST_INFECTION/GlobalDatabase/share/software/Miniconda3/envs.multi-user/vcontact2/bin/mcl"
MCXDUMP="/hwfssz5/ST_INFECTION/GlobalDatabase/share/software/Miniconda3/envs.multi-user/vcontact2/bin/mcxdump"


workingdir="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/test_rec_GPD/testmorehosts" #should be the same place running the script
INPUTdir="/ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/liqian6/GPD/search_by_host"
hostlist="/ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/liqian6/GPD/MCL_test/hostlist.txt"

echo -e "start to prepare fasta files at $(date)"
while read host;
do 

awk -v f="$host" '{OFS="|"; if($1 ~ /^>/) {print $1,f} else {print $0} }' $INPUTdir/${host// /}/${host// /}.phage_sequences.fa >> $workingdir/tmpin.fa
echo -e "done with host $host at $(date)"

done < $hostlist

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

echo -e "start to load abc file at $(date)"

$MCXLOAD  -abc $workingdir/tmpin.abc --stream-mirror  --stream-neg-log10 -stream-tf 'ceil(200)' -write-tab $workingdir/tmpin.tab -o $workingdir/tmpin.mci

echo -e "start to generate clusters at $(date)"

$MCL $workingdir/tmpin.mci -I 1.4

echo -e "start to make tabular results at $(date)"

$MCXDUMP -icl $workingdir/out.tmpin.mci.I14 -tabr  $workingdir/tmpin.tab -o  $workingdir/clusters.mci.I14

echo -e "MCL completed at $(date)"

echo -e "start to extract clusters at $(date)"

mkdir $workingdir/cl

cat $workingdir/clusters.mci.I14 while read line; do echo -e "$line" | tr "\t" "\n" | tr "|" "\t" > $workingdir/tmp ; cut -f1 $workingdir/tmp | sort | uniq | while read contig; do echo -e "$contig|$(awk -v c="$contig" '$1 == c' $workingdir/tmp | cut -f2 | tr "\n" ";")" ; done > $workingdir/cl/cluster_$n.txt ; n=$((n+1))  ; done

rm $workingdir/tmp
echo -e "job completed at $(date)"
