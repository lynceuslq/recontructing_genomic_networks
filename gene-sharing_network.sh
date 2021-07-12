#!/bin/bash

###################################################please define the arguements below######################################################
export PATH="/PATH/TO/PYTHON/:$PATH"
export PATH="/PATH/TO/vcontact2/bin/:$PATH"
CLUSTERONE="/PATH/TO/cluster_one-1.0.jar"
GPDdir="/ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/liqian6/GPD/gut_phage_database"
VCfile="vcset2"#give a name to your output directory
Outfir="/PATH/TO/OUTPUT/${VCfile// /}_vcontact2"
VClist="/PATH/TO/vclist.txt"#a list of viral clusters from GPD to infer gene-sharing networks, should all in intergers and separated by newline elements


##################################################you do not need to change anything below##################################################
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

vcontact2 --raw-proteins $Outfir/tempdir/${VCfile// /}.proteome.faa   --rel-mode 'Diamond' --proteins-fp $Outfir/tempdir/${VCfile// /}.g2g.csv  --db 'ProkaryoticViralRefSeq94-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin $CLUSTERONE --output-dir  $Outfir  

echo -e "finishing vCONTACT2 on $VCfile at $(date)"

grep "vig_" $Outfir/c1.ntw > $Outfir/gpd.c1.ntw

grep "vig_" $Outfir/c1.clusters > $Outfir/gpd.c1.clusters

echo -e "$(cat $Outfir/gpd.c1.clusters | wc -l) gene-sharing clusters from input GPD genomes were found"

n=0
mkdir $Outfir/cl

cut -d "," -f8 $Outfir/gpd.c1.clusters | sed 's/"//g' | while read line; do n=$((n+1)) ; echo -e "cluster${n// /}" >> $Outfir/cl/cl.list ; echo -e "$line" | tr " " "\n" > $Outfir/cl/cluster${n// /}.list ; done

cat $VClist | while read vc ; do cat $Outfir/cl/cl.list | while read cl ; do m=$(grep "VC" $cl.list | wc -l) ; n=$(grep -w "VC${vc}" $Outfir/cl/$cl.list| wc -l) ; echo -e "VC$vc" > $Outfir/cl/VC.$vc.tmp; echo -e "\t$n/$m" >> $Outfir/cl/VC.$vc.tmp; done ; done

echo -e "job completed on $VCfile at $(date)"
