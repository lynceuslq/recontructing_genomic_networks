#!/bin/bash

gpd_cluster="/vcontact2_results/gpd.c1.mo2.clusters"
outpath="/PATH/TO/OUTPUT"
vc2profile="/vcontact2_results/gpd.vConTACT_profiles.csv"
vc2proteins="/vcontact2_results/gpd.vConTACT_proteins.csv"

protsdb="/vcontact2_results/tempdir"
outdir="/PATH/TO/OUTPUTS"
cllist="$outpath/vclist.tmp"
protlist="$outpath/pcvclist.info"
vConTACT_proteins="/vcontact2_results/gpd.vConTACT_proteins.csv"
cdhitprot="/cd-hit-4.8.1/bin/cd-hit"
PLOTCON="/EMBOSS-6.6.0/emboss/plotcon"
MAFFT="mafft/bin/mafft"
RSCRIPT="/PATH/TO/Rscript"
GETDIST="/PATH/TO/get_dist.R"

echo -e "start to process vcontact2 results on $(date)"
n=0

cut -d "," -f8 $gpd_cluster | sed 's/"//g' | while read line; do n=$((n+1)) ; echo -e "$line" | tr " " "\n" > $outpath/cluster${n// /}.list ; echo -e "done with cluster $n at $(date)" ; done

ls $outpath | rev | cut -d "/" -f1 | rev | grep "cluster" | cut -d "." -f1 > $outpath/cl.list

cllist=$outpath/cl.list


cat $cllist | while read cl ; do num=$(echo -e "$cl:$(grep "vig_" $outpath/$cl.list | wc -l )") ; awk -v n="$num" '{OFS="|"; {print $0, n}}' $outpath/$cl.list > $outpath/matchgpdvcon.$cl.list ; echo -e "completed counting on $cl genomes at $(date)"; done

cat $outpath/matchgpdvcon.cluster*.list > $outpath/phage.list

cat $outpath/phage.list | while read ph; do gn=$(echo -e "$ph" | cut -d "|" -f1,2,3) ; echo -e "$ph\t$(grep "$gn" -w $vc2profile | cut -d "," -f2 | sort | tr "\n" ";")" ; done > $outpath/matchvconphageandpc.txt 

cat $vc2proteins | cut -d "," -f4 | tail -n +2 | sort | uniq -c | sort -n -r -k1,1 |sed -e 's/^ *//g' > $outpath/PC.count

cat $outpath/PC.count | grep "PC" | cut -d " " -f2 > $outpath/pc.list

cat $outpath/pc.list | while read pc ; do echo -e "$pc\t$(grep -w "$pc" $outpath/matchvconphageandpc.txt | cut -f1 | cut -d "|" -f4 | sort | uniq -c | sed 's/^ *//g' | sort -k1,1 -n -r | tr "\n" "\t")" ; done > $outpath/matchpctovconvc.nospace.txt

cat $outpath/matchpctovconvc.nospace.txt | cut  -f1 | while read pc; do echo -e "$pc\t$(awk -v p="$pc" -F " " '$2 == p'  $outpath/PC.count | cut -d " " -f1)\t$(awk  -v p="$pc" -F "\t" '$1 == p' $outpath/matchpctovconvc.nospace.txt| cut -f 2-)" ; done > $outpath/matchpctovconvc.nospace.withpccounts.txt 

cat $outpath/matchpctovconvc.nospace.withpccounts.txt  | while read line; do pc=$(echo -e "$line" | cut -f1,2); echo -e "$line" | cut -f 3- | tr "\t" "\n" | while read cl; do echo -e "$pc\t$cl" ; done; done | tr ":" " " | tr " " "\t" > $outpath/matchpctovconvc.pairs.txt 

echo -e "PC\tprot_num_PC\tprot_num_VC\tVC\tgm_num_VC\tprop_PC\tprop_VC" > $outpath/matchpctovconvc.pairs.prop.txt ; cat $outpath/matchpctovconvc.pairs.txt  | awk '{OFS="\t"; if($3 != ""){print $0, $3 / $2, $3 / $5}}' >> $outpath/matchpctovconvc.pairs.prop.txt

cat $outpath/matchpctovconvc.pairs.prop.txt | tail -n +2 | awk '{OFS="\t"; {print $0, $3 / ($2 + $5 - $3)}}' >> $outpath/matchpctovconvc.pairs.stats.txt

cat $outpath/matchpctovconvc.pairs.stats.txt | awk '$8 >= 0.7' | cut -f4 | sort | uniq > $outpath/vclist.tmp 


mkdir $outdir

cat $cllist |grep -v "VC" | while read cl

do

echo "start with $cl at $(date)" 

mkdir $outdir/$cl

cat $protlist | awk -v c="$cl" '$4 == c' | cut -f1 | while read pc; do grep -w "$pc" $vConTACT_proteins | cut -d "|" -f1,2 | sed -e 's/|VC/,/' | tr "," "\t" | while read gn gm vc; do grep -w "$gn" -A1 $protsdb/$vc.proteome.faa > $outdir/$cl/tmpfl; cat $outdir/$cl/tmpfl |  awk -v k="$vc" -v p="$pc" -v c="$cl" '{OFS="|"; if($1 ~ /^>/) {print $0, k , p , c} else {print $0}}' ; done >> $outdir/$cl/$pc.$cl.faa ; done 

echo -e "done with $cl at $(date)" 

done


echo -e "start to cluster selected pc at $(date)"

cat $protlist | cut -f1,4 | tail -n +2 > $outdir/pctovc.list

cat $outdir/pctovc.list | while read pc cl ; do $cdhitprot -i $outdir/$cl/$pc.$cl.faa -o $outdir/$cl/$pc.$cl.clustered.faa -c 0.9 -aS 0.7 -p 1 ; echo -e "completed clustering on $pc from $cl at $(date)"; done

echo -e "start to align and compute conservation status on selected PCs at $(date)"

cat $outdir/pctovc.list | while read pc cl ; do $MAFFT $outdir/$cl/$pc.$cl.faa > $outdir/$cl/$pc.$cl.aligned.faa ; $PLOTCON  -sequences  $outdir/$cl/$pc.$cl.aligned.faa  -winsize 2 -graph  data -goutfile $outdir/$cl/$pc.$cl.aligned.ws2 ;  $PLOTCON  -sequences  $outdir/$cl/$pc.$cl.aligned.faa  -winsize 4  -graph  data -goutfile $outdir/$cl/$pc.$cl.aligned.ws4 ; echo -e "completed clustering on $pc from $cl at $(date)" ; done 

echo -e "start to generate plotcon summary on $(date)"

cat $outdir/pctovc.list | while read pc cl ; do echo -e "$pc\t$cl\t$(cat $outdir/$cl/$pc.$cl.aligned.ws41.dat | grep -v "#" | cut -f2 | tr "\n" ",")" ; done > $outdir/ws4.summary

cat $outdir/pctovc.list | while read pc cl ; do echo -e "$pc\t$cl\t$(cat $outdir/$cl/$pc.$cl.aligned.ws21.dat | grep -v "#" | cut -f2 | tr "\n" ",")" ; done > $outdir/ws2.summary

echo -e "start to generate alignment summary on $(date)"

cp $GETDIST $outdir/get_dist.tmp.R 

sed -i "s/TESTPCVCLIST/$(echo -e "$outdir/pctovc.list" | tr "/" "@" )/" $outdir/get_dist.tmp.R 

sed -i "s/TESTPUTDIR/$(echo -e "$outdir"  | tr "/" "@" )/g" $outdir/get_dist.tmp.R 

cat  $outdir/get_dist.tmp.R  | tr "@" "/" >  $outdir/get_dist.work.R 

rm $outdir/get_dist.tmp.R

$RSCRIPT  $outdir/get_dist.work.R


cat $outdir/pctovc.list |  while read pc vc; do echo -e "$vc\t$pc\t$(head -1 $outdir/$vc/$pc.$vc.dist | tr "\t" "\n" | wc -l)\t$(cat $outdir/$vc/$pc.$vc.dist | tr "\t" "\n" |  grep -v "vig_" | sort -n | awk '{SUM+=$1; nums[NR] = $1}END{print SUM, NR, SUM / NR, nums[1], nums[NR]}')\t$(cat $outdir/$vc/$pc.$vc.clustered.faa  | grep -c ">")" ; done | tr " " "\t" > $outdir/mean_intra_dist.txt 

cat $outdir/pctovc.list |  while read pc vc; do echo -e "$vc\t$pc\t$(cat $outdir/$cl/$pc.$cl.aligned.ws21.dat | grep -v "#" | cut -f2 | wc -l)\t$(cat $outdir/$cl/$pc.$cl.aligned.ws21.dat | grep -v "#" | cut -f2 | awk '$1 >= 2' | wc -l)" ; done $outdir/windowsum.tmp

cat $outdir/windowsum.tmp | awk '{print $0, $4 / $3}' > $outdir/window2.summary

echo -e "job completed at $(date)"
