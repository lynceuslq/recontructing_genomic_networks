#!/bin/bash

#######################################################please define the arguements below######################################################
export PATH="/PATH/TO/ncbi-blast-2.10.1+/bin:$PATH"

MCXLOAD="/PATH/TO/mcxload"
MCL="/PATH/TO/mcl"
MCXDUMP="/PATH/TO/mcxdump"


workingdir="/PATH/TO/WORKING/DIRECTORY/" #should be the same place running the script
INPUTdir="/PATH/TO/GPD/search_by_host" #a directory storing metadata.tsv and phage genomes (phage_sequences.fa) and separated by host species
hostlist="/PATH/TO/hostlist.txt" #a list of host species, should be no space in the names, separated by newline elements


#####################################################you do not need to change anything below###################################################
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

#############################################################making genome clusters##############################################################
echo -e "start to load abc file at $(date)"

$MCXLOAD  -abc $workingdir/tmpin.abc --stream-mirror  --stream-neg-log10 -stream-tf 'ceil(200)' -write-tab $workingdir/tmpin.tab -o $workingdir/tmpin.mci

echo -e "start to generate clusters at $(date)"

$MCL $workingdir/tmpin.mci -I 1.4

echo -e "start to make tabular results at $(date)"

$MCXDUMP -icl $workingdir/out.tmpin.mci.I14 -tabr  $workingdir/tmpin.tab -o  $workingdir/clusters.mci.I14

echo -e "MCL completed at $(date)"


######################################################sorting and remaking genome clusters######################################################

echo -e "start to extract clusters at $(date)"

mkdir $workingdir/cl

cat $workingdir/clusters.mci.I14 while read line; do echo -e "$line" | tr "\t" "\n" | tr "|" "\t" > $workingdir/tmp ; cut -f1 $workingdir/tmp | sort | uniq | while read contig; do echo -e "$contig|$(awk -v c="$contig" '$1 == c' $workingdir/tmp | cut -f2 | tr "\n" ";")" ; done > $workingdir/cl/cluster_$n.txt ; n=$((n+1))  ; done

rm $workingdir/tmp

cut -f1,2,5 $workingdir/tmpin.outfmt6 | tr "|" "\t" | cut -f1,3,5 | sort -u > $workingdir/tmpin.blastn.length.abc

cut -f1,2 $workingdir/tmpin.blastn.length.abc | sort -u | while read c1 c2 ; do awk -v a="$c1" -v b="$c2" '$1 == a && $2 ==b' $workingdir/tmpin.blastn.length.abc | sort -k3,3 -n -r | head -1 ; done > $workingdir/tmpin.blastn.length.sorted.top.abc

echo -e "job completed at $(date)"
