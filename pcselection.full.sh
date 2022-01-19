#!/bin/bash

export PATH="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/miniconda/envs/vContact2/bin:$PATH"
GPDdir="/ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/liqian6/GPD/gut_phage_database"
VCfile="vcset_test"
#change the ouput directory to a writable path
Outfir="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/genesharing/${VCfile// /}_vcontact2"
VClist="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/genesharing/vclist_test.txt"
MAX_OVERLAP="0.2"
DIAMOND="/ldfssz1/ST_META/share/User/zhujie/.conda/envs/bioenv/bin/diamond"


mkdir $Outfir
mkdir $Outfir/tempdir

while read vc; 
do

echo -e "start to files for $vc at $(date)"
awk -v v="$vc" '$3 == v' $GPDdir/GPD_metadata.tsv > $Outfir/tempdir/$vc.metadata.tsv

echo -e "start to extract proteomes for $vc at $(date)"
cat $Outfir/tempdir/$vc.metadata.tsv | cut -f1 | while read con; do grep ">${con// /}_" -A1 $GPDdir/GPD_proteome.faa ; done > $Outfir/tempdir/$vc.proteome.faa

grep  ">" $Outfir/tempdir/$vc.proteome.faa | sed -e 's/>//' | while read id ; do genome=$(echo -e "$id" | cut -d "_" -f1,2) ; ph=$(grep -w "$genome" $Outfir/tempdir/$vc.metadata.tsv | cut -f15) ; echo -e "$id,$genome|VC${vc}|$ph,VC${vc}"; done > $Outfir/tempdir/$vc.g2g.csv


echo -e "prepared file for $vc at $(date)"

done < $VClist

echo -e "protein_id,contig_id,keywords" > $Outfir/tempdir/${VCfile// /}.g2g.csv
cat $VClist | while read vc; do cat $Outfir/tempdir/$vc.g2g.csv >> $Outfir/tempdir/${VCfile// /}.g2g.csv ; done

cat $VClist | while read vc; do cat $Outfir/tempdir/$vc.proteome.faa >> $Outfir/tempdir/${VCfile// /}.proteome.faa ; done

echo -e "starting to contsruct network at $(date)"

vcontact2 --raw-proteins $Outfir/tempdir/${VCfile// /}.proteome.faa   --rel-mode 'Diamond' --proteins-fp $Outfir/tempdir/${VCfile// /}.g2g.csv  --db 'ProkaryoticViralRefSeq94-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin /jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/tools/cluster_one-1.0.jar --output-dir  $Outfir  

#here viral clusters are reconstructed with clusterone with maximum overlapping rate set to 0.2
java -jar /jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/tools/cluster_one-1.0.jar  $Outfir/c1.ntw  --input-format edge_list --output-format csv --min-density 0.3 --min-size 2 --max-overlap $MAX_OVERLAP --penalty 2.0 --haircut 0.55 --merge-method single --similarity match --seed-method nodes > $Outfir/c1.mo2.clusters

echo -e "finishing vCONTACT2 on $VCfile at $(date)"

grep "vig_" $Outfir/c1.ntw > $Outfir/gpd.c1.ntw

grep "vig_" $Outfir/c1.clusters > $Outfir/gpd.c1.clusters

grep "vig_" $Outfir/c1.mo2.clusters > $Outfir/gpd.c1.mo2.clusters

grep "vig_" $Outfir/vConTACT_profiles.csv > $Outfir/gpd.vConTACT_profiles.csv

grep "vig_" $Outfir/vConTACT_proteins.csv > $Outfir/gpd.vConTACT_proteins.csv

echo -e "$(cat $Outfir/gpd.c1.mo2.clusters | wc -l) gene-sharing clusters from input GPD genomes were found"

echo -e "VC reconstruction completed for $VCfile at $(date)"


gpd_cluster=$Outfir/gpd.c1.mo2.clusters
outpath=$Outfir/cl2
vc2profile=$Outfir/gpd.vConTACT_profiles.csv
vc2proteins=$Outfir/gpd.vConTACT_proteins.csv

mkdir $outpath

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


#computing jaccard index and select pc passing the threshold
cat $outpath/matchpctovconvc.pairs.prop.txt | tail -n +2 | awk '{OFS="\t"; {print $0, $3 / ($2 + $5 - $3)}}' >> $outpath/matchpctovconvc.pairs.stats.txt

#jaccard index as 0.7
cat $outpath/matchpctovconvc.pairs.stats.txt | awk '$8 >= 0.7' | cut -f4 | sort | uniq > $outpath/vclist.tmp 

echo -e "PC\tprot_num_PC\tprot_num_VC\tVC\tgm_num_VC\tprop_PC\tprop_VC\tjaccard" > $outpath/pcvclist.info 

#jaccard index as 0.7
awk '$8 >= 0.7' $outpath/matchpctovconvc.pairs.stats.txt  >> $outpath/pcvclist.info

echo -e "splitting viral clusters completed at $(date)"


#defining input and output directories
protsdb=$Outfir/tempdir
outdir=$outpath/selected_pc
cllist=$outpath/vclist.tmp
protlist=$outpath/pcvclist.info
vConTACT_proteins=$Outfir/gpd.vConTACT_proteins.csv


#defining the path of softwares here
cdhitprot="/ldfssz1/ST_META/share/User/zhujie/.conda/envs/bioenv/bin/cd-hit"
PLOTCON="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/tools/EMBOSS-6.6.0/emboss/plotcon"
MAFFT="/ldfssz1/ST_META/share/User/zhujie/.conda/envs/bioenv/bin/mafft"
RSCRIPT="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/miniconda/bin/Rscript"

mkdir $outdir

cat $cllist |grep -v "VC" | while read cl

do

echo "start with $cl at $(date)" 

mkdir $outdir/$cl

#loading protein metadata and sequences
cat $protlist | awk -v c="$cl" '$4 == c' | cut -f1 | while read pc; do grep -w "$pc" $vConTACT_proteins | cut -d "|" -f1,2 | sed -e 's/|VC/,/' | tr "," "\t" | while read gn gm vc; do grep -w "$gn" -A1 $protsdb/$vc.proteome.faa > $outdir/$cl/tmpfl; cat $outdir/$cl/tmpfl |  awk -v k="$vc" -v p="$pc" -v c="$cl" '{OFS="|"; if($1 ~ /^>/) {print $0, k , p , c} else {print $0}}' ; done >> $outdir/$cl/$pc.$cl.faa ; done 

echo -e "done with $cl at $(date)" 

done


echo -e "start to cluster selected pc at $(date)"

cat $protlist | cut -f1,4 | tail -n +2 > $outdir/pctovc.list

#performing clutstering
cat $outdir/pctovc.list | while read pc cl ; do $cdhitprot -i $outdir/$cl/$pc.$cl.faa -o $outdir/$cl/$pc.$cl.clustered.faa -c 0.9 -aS 0.7 -p 1 ; echo -e "completed clustering on $pc from $cl at $(date)"; done

echo -e "start to align and compute conservation status on selected PCs at $(date)"

#generating alignments and calculate conservation status by wondow sliding
cat $outdir/pctovc.list | while read pc cl ; do $MAFFT $outdir/$cl/$pc.$cl.faa > $outdir/$cl/$pc.$cl.aligned.faa ; $PLOTCON  -sequences  $outdir/$cl/$pc.$cl.aligned.faa  -winsize 2 -graph  data -goutfile $outdir/$cl/$pc.$cl.aligned.ws2 ;  $PLOTCON  -sequences  $outdir/$cl/$pc.$cl.aligned.faa  -winsize 4  -graph  data -goutfile $outdir/$cl/$pc.$cl.aligned.ws4 ; echo -e "completed clustering on $pc from $cl at $(date)" ; done 

echo -e "start to generate plotcon summary on $(date)"

#summarizing conservation status by site
cat $outdir/pctovc.list | while read pc cl ; do echo -e "$pc\t$cl\t$(cat $outdir/$cl/$pc.$cl.aligned.ws41.dat | grep -v "#" | cut -f2 | tr "\n" ",")" ; done > $outdir/ws4.summary

cat $outdir/pctovc.list | while read pc cl ; do echo -e "$pc\t$cl\t$(cat $outdir/$cl/$pc.$cl.aligned.ws21.dat | grep -v "#" | cut -f2 | tr "\n" ",")" ; done > $outdir/ws2.summary

echo -e "start to generate alignment summary on $(date)"

cp /jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/genesharing/vcset6_split_vcontact2/get_dist.R $outdir/get_dist.tmp.R 

sed -i "s/TESTPCVCLIST/$(echo -e "$outdir/pctovc.list" | tr "/" "@" )/" $outdir/get_dist.tmp.R 

sed -i "s/TESTPUTDIR/$(echo -e "$outdir"  | tr "/" "@" )/g" $outdir/get_dist.tmp.R 

cat  $outdir/get_dist.tmp.R  | tr "@" "/" >  $outdir/get_dist.work.R 

rm $outdir/get_dist.tmp.R

#calculating pairwise distance
$RSCRIPT  $outdir/get_dist.work.R


cat $outdir/pctovc.list |  while read pc vc; do echo -e "$vc\t$pc\t$(head -1 $outdir/$vc/$pc.$vc.dist | tr "\t" "\n" | wc -l)\t$(cat $outdir/$vc/$pc.$vc.dist | tr "\t" "\n" |  grep -v "vig_" | sort -n | awk '{SUM+=$1; nums[NR] = $1}END{print SUM, NR, SUM / NR, nums[1], nums[NR]}')\t$(cat $outdir/$vc/$pc.$vc.clustered.faa  | grep -c ">")" ; done | tr " " "\t" > $outdir/mean_intra_dist.txt 

cat $outdir/pctovc.list |  while read pc cl; do echo -e "$cl\t$pc\t$(cat $outdir/$cl/$pc.$cl.aligned.ws21.dat | grep -v "#" | cut -f2 | wc -l)\t$(cat $outdir/$cl/$pc.$cl.aligned.ws21.dat | grep -v "#" | cut -f2 | awk '$1 >= 2' | wc -l)" ; done >  $outdir/windowsum.tmp

cat $outdir/windowsum.tmp | awk '{print $0, $4 / $3}' | tr " " "\t" > $outdir/windows.selected.list

#selecting pc by sequential distance
cat $outdir/mean_intra_dist.txt  | awk '$6 <= 0.25 && $8 <= 0.8 && $ 9 <= 2' | cut -f1,2 > $outdir/tmplist2
#selecting pc by conservation plots
cat  $outdir/windows.selected.list | awk '$5 >= 0.85' |cut -f1,2 >  $outdir/tmplist1

cat $outdir/tmplist2 $outdir/tmplist1  | sort | uniq -d > $outdir/selected_vcpc.list

#loading and reformatting selected sequences
cat $outdir/selected_vcpc.list | while read vc pc ; do echo -e "$(cat $outdir/$vc/$pc.$vc.clustered.faa   | tr "\n" "\t")" |   tr ">" "\n" | sed '/^[[:space:]]*$/d' | while read c1 c2; do a=$(echo -e "$c2" | wc -c) ; b=$(($a-1)) ; echo -e ">$c1|${b}aa\n$c2" ; done ; done >  $outdir/selected_pc.len.clstr.faa 

cat $outdir/selected_vcpc.list  | cut -f1 | sort | uniq | while read vc ; do cat $outpath/phage.list | grep -w "$vc" | grep "vig_" >> $outdir/selected_pc.phage.list ; echo -e "done with $vc" ; done



preselection="$outdir/preselectionprots.2.list"

#setting the minimum amount of PCs for a VC
cat $outdir/selected_vcpc.list | cut -f1 | uniq -c | awk '$1 > 8' | rev | cut -d " " -f1 | rev > $outdir/selected.r1.list

grep -v -w -f $outdir/selected.r1.list $outpath/pcvclist.info | awk '{OFS="\t"; {print $4, $1}}' > $outdir/tmplist1

cat $outdir/selected_vcpc.list | cut -f2 > $outdir/tmplist2

grep -v -w -f $outdir/tmplist2 $outdir/tmplist1 > $preselection



export PATH="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/miniconda/bin:$PATH"

indir=$outdir
workingdir=$outdir
outdir="$indir/preselection"

tail -n +2 $outpath/pcvclist.info | cut -f4 | sort | uniq > $workingdir/allcl2sig.list 

mkdir $outdir
echo -e "job started on $(date)"

cat $preselection | cut -f1,2 | while read vc pc
do

cp $indir/$vc/$pc.$vc.aligned.faa $outdir/$pc.$vc.aligned.faa

Gblocks $outdir/$pc.$vc.aligned.faa  -t=p -e=-gb1 -b4=20

$cdhitprot -i $outdir/$pc.$vc.aligned.faa-gb1  -o $outdir/$pc.$vc.aligned.faa-gb1.clustered.faa -c 0.8 -aS 0.7 -p 1

echo -e "completed with $vc $pc"

done

echo -e "job completed at $(date)"

grep -c ">" $outdir/*.aligned.faa-gb1.clustered.faa | tr ":" "\t" | awk '$2 > 0 && $2 <= 3' | cut -f1 | while read file ; do cat $file ; done > $workingdir/rtwoselected.clustered.faa

cat $workingdir/rtwoselected.clustered.faa  | tr "\n" "\t" | tr ">" "\n" | sed '/^[[:space:]]*$/d' | while read line; do echo -e ">$(echo -e "$line" | cut -f1)\n$(echo -e "$line" | cut -f2- | sed -e "s/\t//g")" ; done > $workingdir/moreselected.clustered.faa


echo -e ">" >> $workingdir/moreselected.clustered.faa
cat $workingdir/moreselected.clustered.faa  | tr "\n" "\t" | tr ">" "\n" | sed '/^[[:space:]]*$/d' | while read acc seq ; do newacc=$(echo -e "$acc" | cut -d "|" -f1,2,3,4) ; a=$(echo -e "$seq" | wc -c) ; b=$((a-1)) ;  echo -e ">$newacc|${b}aa\n$seq" ; done  > $workingdir/moreselected.clustered.len.faa

grep ">"  $workingdir/moreselected.clustered.len.faa | tr "|" "\t" | awk '$5 >= 15' | cut -f3,4 |  sort | uniq > $workingdir/newtoadd.pcvc.list


cat $workingdir/newtoadd.pcvc.list |  while read pc vc ; do cat $outdir/$pc.$vc.aligned.faa-gb1.clustered.faa ; done > $workingdir/newtoadd.filtered.clustered.faa


cat $workingdir/newtoadd.filtered.clustered.faa | tr "\n" "\t" | tr ">" "\n" | sed '/^[[:space:]]*$/d' | while read line; do echo -e ">$(echo -e "$line" | cut -f1)\n$(echo -e "$line" | cut -f2- | sed -e "s/\t//g")" ; done | sed -e "s/ //g" > $workingdir/newtoadd.filtered.clustered.formatted.faa

cat $workingdir/selected_pc.len.clstr.faa $workingdir/newtoadd.filtered.clustered.formatted.faa > $Outfir/total.selected_pc.len.clstr.faa

echo -e "start to remove inter-homo PCs at $(date)"

$DIAMOND makedb -d $Outfir/total.selected_pc.len.clstr --in $Outfir/total.selected_pc.len.clstr.faa

$DIAMOND blastx -d $Outfir/total.selected_pc.len.clstr -q $GPDdir/GPD_sequences.fa -f 6 -o $Outfir/pairwiseblast.tab --threads 2

cat $Outfir/pairwiseblast.tab |  awk -F "\t" '$3 >= 75 && $4 >= 15 ' | cut -f1,2 | tr "|" "\t" | cut -f1,4,5 | sort | uniq | while read gm pc cl ; do echo -e "$gm\t$(grep -w "$gm" $outpath/phage.list | tr ":" "|" | cut -d "|" -f4 | tr "\n" ":")\t$pc\t$cl"  ; done | awk -F "\t" '$2 != ""' | while read gm cl1 pc cl2 ; do if [[ "$cl1" != *"$cl2:"* ]]; then echo -e "$pc"; fi ; done | sort | uniq > $Outfir/pairwiseblast.interhomo.tab

cat $Outfir/total.selected_pc.len.clstr.faa | grep ">" | tr "|" "\t" | cut -f3,4 | sort | uniq > $Outfir/preremoval.pcvc.list

grep -v -w -f $Outfir/pairwiseblast.interhomo.tab $Outfir/preremoval.pcvc.list > $Outfir/final.pcvc.list

cat $Outfir/final.pcvc.list  | while read pc vc ; do grep -w "$pc" -A1 $Outfir/total.selected_pc.len.clstr.faa ; done > $Outfir/final.selected_pc.len.clstr.faa

grep ">" $Outfir/final.selected_pc.len.clstr.faa | cut -d "|" -f4 | sort | uniq | while read vc ; do cat $outpath/phage.list | grep -w "$vc" | grep "vig_" >> $Outfir/final.selected_pc.phage.list ; echo -e "done with $vc" ; done


echo -e "job completed at $(date)"
