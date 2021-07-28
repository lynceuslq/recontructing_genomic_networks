#!/bin/bash

gpd_cluster="/PATH/TO/VCONTACT2/RESULTS/c1.clusters"
outpath="/PATH/TO/VCONTACT2/CLUSTERS"
cllist="/PATH/TO/VCONTACT2/CLUSTERS/cl.list"
vc2profile="/PATH/TO/VCONTACT2/RESULTS/gpd.vConTACT_profiles.csv"
vc2proteins="/PATH/TO/VCONTACT2/RESULTS/gpd.vConTACT_proteins.csv"

echo -e "start to process vcontact2 results on $(date)"
n=0

cut -d "," -f8 $gpd_cluster | sed 's/"//g' | while read line; do n=$((n+1)) ; echo -e "$line" | tr " " "\n" > $outpath/cluster${n// /}.list ; echo -e "done with cluster $n at $(date)" ; done

cat $cllist | while read cl ; do num=$(echo -e "$cl:$(grep "vig_" $outpath/$cl.list | wc -l )") ; awk -v n="$num" '{OFS="|"; {print $0, n}}' $outpath/$cl.list > $outpath/matchgpdvcon.$cl.list ; echo -e "completed counting on $cl genomes at $(date)"; done

cat $outpath/matchgpdvcon.cluster*.list > $outpath/phage.list

cat $outpath/phage.list | while read ph; do gn=$(echo -e "$ph" | cut -d "|" -f1,2,3) ; echo -e "$ph\t$(grep "$gn" -w $vc2profile | cut -d "," -f2 | sort | tr "\n" ";")" ; done > $outpath/matchvconphageandpc.txt 

cat $vc2proteins | cut -d "," -f4 | tail -n +2 | sort | uniq -c | sort -n -r -k1,1 |sed -e 's/^ *//g' > $outpath/PC.count

cat $outpath/pc.list | while read pc ; do echo -e "$pc\t$(grep -w "$pc" $outpath/matchvconphageandpc.txt | cut -f1 | cut -d "|" -f4 | sort | uniq -c | sed 's/^ *//g' | sort -k1,1 -n -r | tr "\n" "\t")" ; done > $outpath/matchpctovconvc.nospace.txt

cat $outpath/matchpctovconvc.nospace.txt | cut  -f1 | while read pc; do echo -e "$pc\t$(awk -v p="$pc" -F " " '$2 == p'  $outpath/PC.count | cut -d " " -f1)\t$(awk  -v p="$pc" -F "\t" '$1 == p' $outpath/matchpctovconvc.nospace.txt| cut -f 2-)" ; done > $outpath/matchpctovconvc.nospace.withpccounts.txt 

cat $outpath/matchpctovconvc.nospace.withpccounts.txt  | while read line; do pc=$(echo -e "$line" | cut -f1,2); echo -e "$line" | cut -f 3- | tr "\t" "\n" | while read cl; do echo -e "$pc\t$cl" ; done; done | tr ":" " " | tr " " "\t" > $outpath/matchpctovconvc.pairs.txt 

echo -e "PC\tprot_num_PC\tprot_num_VC\tVC\tgm_num_VC\tprop_PC\tprop_VC" > $outpath/matchpctovconvc.pairs.prop.txt ; cat $outpath/matchpctovconvc.pairs.txt  | awk '{OFS="\t"; if($3 != ""){print $0, $3 / $2, $3 / $5}}' >> $outpath/matchpctovconvc.pairs.prop.txt


echo -e "job completed at $(date)"
