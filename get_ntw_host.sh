#!/bin/bash

##############################################please deifne the arguements below################################################

BLASTN_tab="/PATH/TO/blastn.length.sorted.top.abc"
#can be found in the outputs of blastnandcluster.v2.sh 

PHAGEANDHOST="/PATH/TO/genomewithhost.list"
#a file of phage genomes and their host species seperated by "|", if there are multiple hosts for the same phage, hosts should be separated by ";"
#example ivig_1019|Clostridium_M_bolteae;Dorea_longicatena;

OUTPUTPATH="/PATH/TO/OUTPUT/culsterbyhost"
hostlist="/PATH/TO/hostlist.txt"


##############################################you do not need to change anything below##########################################
cat $hostlist | while read host
do

echo -e "started with host $host at $(date)"

grep "$host" $PHAGEANDHOST | cut -d "|" -f1 | while read con; do grep "$con" -w $BLASTN_tab ; done > $OUTPUTPATH/$host.test.list

cat $OUTPUTPATH/$host.test.list |  awk '{print log($3)}' > $OUTPUTPATH/tmpc1

cut -f1 $OUTPUTPATH/$host.test.list | while read con; do awk -F "|" -v c="$con" '$1 == c' $PHAGEANDHOST | cut -d "|" -f2 ; done > $OUTPUTPATH/tmpc2

cut -f2 $OUTPUTPATH/$host.test.list | while read con; do awk -F "|" -v c="$con" '$1 == c' $PHAGEANDHOST | cut -d "|" -f2 ; done > $OUTPUTPATH/tmpc3

paste $OUTPUTPATH/$host.test.list $OUTPUTPATH/tmpc1 $OUTPUTPATH/tmpc2 $OUTPUTPATH/tmpc3 > $OUTPUTPATH/$host.gpd.ntw

echo -e "completetd host $host at $(date)"

done
echo -e "job completed at $(date)"
