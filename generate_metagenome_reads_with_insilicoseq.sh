#!/bin/bash

BEDTOOLS="/PATH/TO/bedtools2/bin/bedtools"
bacgmdir="/PATH/TO/bacterial_hosts_genomes"
gpdfna="/PATH/TO/GPD/gut_phage_database/GPD_sequences.fa"
gpdmetadata="/PATH/TO/GPD/gut_phage_database/GPD_metadata.tsv"
phage2hosts="/PATH/TO/matchphageandisolates.txt"
hostmetadata="/PATH/TO/bacterialhost.concise.tsv"


helpFunction()
{
   echo ""
   echo "Usage: $0 -d workingdir -s metagenome_names -n number_of_metagenomes -h to_add_bacterial_hosts_or_not -m model_of_simulation -a number_of_reads"
   echo "header of tabular outputs: qacc sacc qlen slen pident evalue bitscore mismatch staxids sscinames scomnames sblastnames sskingdoms"
   echo -e "\t-d drectory to outputs"
   echo -e "\t-s give a name to your metagenome sets"
   echo -e "\t-n number of metagenomes to simulate"
   echo -e "\t-h true if include genomes of bacterial hosts, false to exclude them"
   echo -e "\t-m model of insilicoseq simulations, eg. novaseq, miseq, hiseq"
   echo -e "\t-a Number of reads to generate, eg. 1000000. Allows suffixes k, K, m, M, g and G (ex 0.5M for 500000)"
   exit 1 # Exit script after printing help
}

while getopts "t:d:q:e:o:" opt
do
   case "$opt" in
      d ) workingdir="$OPTARG" ;;
      s ) metagenome_names="$OPTARG" ;;
      n ) number_of_metagenomes="$OPTARG" ;;
      h ) to_add_bacterial_hosts_or_not="$OPTARG" ;;
      m ) model_of_simulation="$OPTARG" ;;
      a ) number_of_reads="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done


if [ -z "$workingdir" ] || [ -z "$metagenome_names" ] || [ -z "$number_of_metagenomes" ] || [ -z "$to_add_bacterial_hosts_or_not" ] || [ -z "$number_of_reads" ] || [ -z "$model_of_simulation" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
else

seq 1 1 $number_of_metagenomes | awk -v s="$metagenome_names" '{OFS="_"; {print "metaset",s, $1}}' | while read setacc; 
do
out_acc=$setacc

cut -f1,3,15 $gpdmetadata | tr "\t" "|" | tail -n +2 | shuf | head -500 > $workingdir/$out_acc.ph.list

phgmlist=$workingdir/$out_acc.ph.list
echo -e "start to extract phage genome set $out_acc at $(date)"
cat $phgmlist | cut -d "|" -f1 | while read acc ; do grep -w "$acc" -A1 $gpdfna >> $workingdir/$out_acc.fna ; done

#searching for bacterial host for target genomes and adding top 5 most common ones to metagenomes
if [ $to_add_bacterial_hosts_or_not == "true" ]
then
cat $workingdir/$out_acc.ph.list  | cut -d "|" -f1 | while read gm ; do awk -v g="$gm" '$1 == g' $phage2hosts ; done | cut -f2 | tr "," "\n" | sort | uniq -c | sort -n -r -k1,1 | head  -50 | sed -e "s/^ *//g" | cut -d " " -f2 | while read acc; do grep -w "$acc" $hostmetadata ; done  | cut -f3 | head -5  > $workingdir/$out_acc.bac.list

bacgmlist=$workingdir/$out_acc.bac.list

echo -e "start to extract bacterial genomes set $out_acc at $(date)"

cat $bacgmlist | while read gm ; do cat $bacgmdir/$gm.fna >> $workingdir/$out_acc.fna ; done
fi


echo -e "completed extraction at $(date)"

echo -e "Start to generate paired end reads at $(date)"

grep ">" $workingdir/$out_acc.fna | sed -e 's/>//' | cut -f1  | while read gm; do echo -e "$gm\t$(shuf -i 1-100 -n 1)" ; done | awk -F "\t" '{OFS="\t"; {if($1 ~ /^GUT/) {print $1, $2 * 8} else {print $0}}}' > $workingdir/$out_acc.tmp

sum=$(awk '{SUM+=$2}END{print SUM}' $workingdir/$out_acc.tmp)

cat $workingdir/$out_acc.tmp | awk -v s=$sum '{OFS="\t"; {print $1, $2 / s}}' > $workingdir/${out_acc}_abundance.txt

iss generate --genomes $workingdir/$out_acc.fna  --model ${model_of_simulation/ //} --output $workingdir/$out_acc --abundance_file $workingdir/${out_acc}_abundance.txt –n_reads ${number_of_reads/ //}

echo -e "job completed with set $setacc at $(date), outputting metagenome reads $workingdir/${out_acc}_R*.fastq"
done
fi
