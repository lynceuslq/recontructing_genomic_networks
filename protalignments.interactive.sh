#!/bin/bash

######################################################defining arguements here###################################################################
REFORMATTER="/PATH/TO/reformat.sh"
DIAMOND="/PATH/TO/Diamond/2.0.9/diamond"
DIAMOND_DB="/PATH/TO/cl2_recoverprots.clustered.len.overlappingremoved.dmnd"
BEDTOOLS="/PATH/TO/bedtools2/bin/bedtools"
PROTLIST="//PATH/TO/cl2_recoverprots.clustered.len.overlappingremoved.length.txt"
phlist="/PATH/TO/gpdvcmorethan3.clustered.len.overlappingremoved.phage.list"
RSCRIPT="/PATH/TO/Rscript"
SUMUPR="/PATH/TO/sumup.protalign.R"

#####################################################you do not need to change anything below###################################################
helpFunction()
{
   echo ""
   echo "Usage: $0 -1 read1.fastq.gz -2 read1.fastq.gz -p parameters -o output_directory -s sample_accession -t threads"
   echo -e "\t-1 the path to read 1 of a sample, or the path to the fastq file if the input is SE"
   echo -e "\t-2 the path to read 2 of a sample, or leave the option blank if the input fastq is SE"
   echo -e "\t-p parameter selection, please select from sensitive, intermediate and specific"
   echo -e "\t-o the path to the output directory"
   echo -e "\t-s sample accession"
   echo -e "\t-t threads  for diamond blastx"
   exit 1 # Exit script after printing help
}

while getopts "1:2:p:o:s:t:" opt
do
   case "$opt" in
      1 ) READ1="$OPTARG" ;;
      2 ) READ2="$OPTARG" ;;
      p ) kvalue="$OPTARG" ;;
      o ) OUTDIR1="$OPTARG" ;;
      s ) ACCESSION="$OPTARG" ;;
      t ) THREADS="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done


#kvalue="sensitive"
OUTDIR="$OUTDIR1/$kvalue"

if [ -z "$READ1" ] || [ -z "$kvalue" ] || [ -z "$OUTDIR1" ] || [ -z "$ACCESSION" ] || [ -z "$THREADS" ] 
then
   echo "Some or all of the parameters are empty";
   helpFunction

else
mkdir $OUTDIR1
mkdir $OUTDIR
echo -e "start to process fastq reads from $ACCESSION at $(date)" 
mkdir $OUTDIR/$ACCESSION

OUTPATH="$OUTDIR/$ACCESSION"

if [ -z "$READ2" ]
then

zcat $READ1 > $OUTPATH/${ACCESSION// /}.tmp.fq

else

echo -e "start to merge paired end reads of $ACCESSION at $(date)"

$REFORMATTER in1=$READ1  in2=$READ2  out=$OUTPATH/${ACCESSION// /}.tmp.fq

echo -e "fninishing merging paired end reads of $ACCESSION at $(date), results stored as $OUTPATH/$ACCESSION.tmp.fq"
fi


echo -e "start to align reads of $ACCESSION at $(date)"

$DIAMOND blastx  -d  $DIAMOND_DB  -q $OUTPATH/$ACCESSION.tmp.fq -o $OUTPATH/$ACCESSION.tab -f 6 -k 10 --threads $THREADS

echo -e "finishing alignment of $OUTPATH/$ACCESSION.tmp.fq, start to genrate coverage at $(date)"


if [ $kvalue == "sensitive" ]
then
echo -e "using $kvalue parameters at $(date)"
cat $OUTPATH/$ACCESSION.tab | awk '$3 >= 85 || $12 >= 120'  | awk '{if($9 < $10){print($2"\t"$9-1"\t"$10)}else{print($2"\t"$10-1"\t"$9)}}' |  sort -k1,1 -k2,2n > $OUTPATH/$ACCESSION.filtered.bed


elif [ $kvalue == "intermediate" ]
then
echo -e "using $kvalue parameters at $(date)"
cat $OUTPATH/$ACCESSION.tab | awk '$3 >= 85 || $12 >= 150'  | awk '{if($9 < $10){print($2"\t"$9-1"\t"$10)}else{print($2"\t"$10-1"\t"$9)}}' |  sort -k1,1 -k2,2n > $OUTPATH/$ACCESSION.filtered.bed


elif [ $kvalue == "specific" ]
then
echo -e "using $kvalue parameters at $(date)"
cat $OUTPATH/$ACCESSION.tab | awk '$3 >= 90 || $12 >= 200'  | awk '{if($9 < $10){print($2"\t"$9-1"\t"$10)}else{print($2"\t"$10-1"\t"$9)}}' |  sort -k1,1 -k2,2n > $OUTPATH/$ACCESSION.filtered.bed


else
echo -e "exit with invalidate parameter selection"
fi


$BEDTOOLS  genomecov -i $OUTPATH/$ACCESSION.filtered.bed -g $PROTLIST > $OUTPATH/$ACCESSION.filtered.cov.txt 

rm $OUTPATH/$ACCESSION.tmp.fq

echo -e "finishing alignment on $ACCESSION at $(date)"


echo -e "summarizing on $ACCESSION at $(date)"

$RSCRIPT --vanilla $SUMUPR $ACCESSION $kvalue $phlist $OUTDIR1 $OUTPATH

echo -e "summarizing on $ACCESSION completed at $(date)"

fi

echo -e "job completed at $(date)"
