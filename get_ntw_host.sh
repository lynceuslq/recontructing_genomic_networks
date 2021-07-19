#!/bin/bash

##############################################please deifne the arguements below################################################

BLASTN_tab="/PATH/TO/blastn.length.sorted.top.abc"
#can be found in the outputs of blastnandcluster.v2.sh 

GENOMELIST="/PATH/TO/genomewithhost.list"
#a file of phage genomes and their host species seperated by "|", if there are multiple hosts for the same phage, hosts should be separated by ";"
#example ivig_1019|Clostridium_M_bolteae;Dorea_longicatena;

METADATA="/PATH/TO/GPD_metadata.tsv"
HOSTBLAST="/PATH/TO/HOSTblastn.tab"
OUTPUTPATH="/PATH/TO/OUTPUT/culsterbyhost"
hostlist="/PATH/TO/hostlist.txt"


##############################################you do not need to change anything below##########################################
cat $hostlist | while read host hostacc
do

echo -e "started with host $host at $(date)"

grep "$host" $GENOMELIST | cut -f1 | while read con; do grep "$con" -w $BLASTN_tab ; done > $OUTPUTPATH/$host.test.list

cat $OUTPUTPATH/$host.test.list |  awk '{print log($3)}' > $OUTPUTPATH/tmpc1

cut -f1 $OUTPUTPATH/$host.test.list | while read con; do awk  -v c="$con" '$1 == c' $GENOMELIST | cut  -f2 ; done > $OUTPUTPATH/tmpc2

cut -f1 $OUTPUTPATH/$host.test.list | while read con; do awk -F "\t" -v c="$con" '$1 == c' $METADATA | cut -f3,15 ; done > $OUTPUTPATH/tmpc3

cut -f2 $OUTPUTPATH/$host.test.list | while read con; do awk  -v c="$con" '$1 == c' $GENOMELIST | cut  -f2 ; done > $OUTPUTPATH/tmpc4

cut -f2 $OUTPUTPATH/$host.test.list | while read con; do awk -F "\t" -v c="$con" '$1 == c' $METADATA | cut -f3,15 ; done > $OUTPUTPATH/tmpc5

paste $OUTPUTPATH/$host.test.list $OUTPUTPATH/tmpc1 $OUTPUTPATH/tmpc2 $OUTPUTPATH/tmpc3 $OUTPUTPATH/tmpc4 $OUTPUTPATH/tmpc5 > $OUTPUTPATH/$host.gpd.ntw

cat $HOSTBLAST/$hostacc.tmp.out | sort -n -k7,7 | cut -f 1-8 | awk '$5 >= 500' | cut -f1,2 | uniq | while read c1 c2 ; do awk -v a="$c1" -v b="$c2" '$1 == a && $2 ==b' $HOSTBLAST/$hostacc.tmp.out | sort -k5,5 -n -r | head -1 ; done | awk -v h="$host" '{OFS="\t"; {print $1, $2, $5, log($5), h, "host", "host"}}' > $HOSTBLAST/$hostacc.tmp.out.sorted

cut -f2 $HOSTBLAST/$hostacc.tmp.out.sorted | while read con; do awk  -v c="$con" '$1 == c' $GENOMELIST | cut  -f2 ; done > $OUTPUTPATH/tmpc4 
cut -f2 $HOSTBLAST/$hostacc.tmp.out.sorted | while read con; do awk -F "\t" -v c="$con" '$1 == c' $METADATA | cut -f3,15 ; done > $OUTPUTPATH/tmpc5

paste $HOSTBLAST/$hostacc.tmp.out.sorted $OUTPUTPATH/tmpc4  $OUTPUTPATH/tmpc5 >> $OUTPUTPATH/$host.gpd.ntw


echo -e "completetd host $host at $(date)"

done
echo -e "job completed at $(date)"
