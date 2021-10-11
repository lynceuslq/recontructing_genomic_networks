#!/bin/bash

kvalue="adjsen"
INPATH="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/abd_dist_model_testsets/onlyclusteredphage"
OUTDIR1="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/recoverprots_gpdvcmorethan3_filterprotcov"
OUTDIR="$OUTDIR1/$kvalue"
REFORMATTER="/hwfssz5/ST_INFECTION/GlobalDatabase/user/fengqikai/software/bbmap/reformat.sh"
DIAMOND="/hwfssz5/ST_INFECTION/GlobalDatabase/user/fengqikai/software/Diamond/2.0.9/diamond"
DIAMOND_DB="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/recoverprots_gpdvcmorethan3/recoverprots_gpdvcmorethan3.dmnd"
BEDTOOLS="/zfssz2/ST_MCHRI/COHORT/fengqikai/software/bedtools2/bedtools2/bin/bedtools"
PROTLIST="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/recoverprots_gpdvcmorethan3/recoverprots_gpdvcmorethan3.protlength.tab"
ABUNDANCEFILES="$OUTDIR1/testabd"
ACCESSIONLIST="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/recoverprots_gpdvcmorethan3/SAMPLELIST"
phlist="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/recoverprots_gpdvcmorethan3/gpdvcmorethan3_vcontact2.cl2_recoverprots.phage.list"
SIMULATIONDIR="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/abd_dist_model_testsets/abd_pack"
RSCRIPT="/hwfssz5/ST_INFECTION/HPV/ouzhihua/program/miniconda2/bin/Rscript"
SUMUPR="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/recoverprots_gpdvcmorethan3/sumup.argparse.protcovfilter.R"

mkdir $OUTDIR
mkdir $ABUNDANCEFILES

cat $ACCESSIONLIST | while read ACCESSION;
do

echo -e "start to process fastq reads from $ACCESSION at $(date)" 
mkdir $OUTDIR/$ACCESSION

OUTPATH="$OUTDIR/$ACCESSION"

echo -e "start to merge paired end reads of $ACCESSION at $(date)"

$REFORMATTER in1=$INPATH/${ACCESSION// /}_R1.fastq  in2=$INPATH/${ACCESSION// /}_R2.fastq  out=$OUTPATH/${ACCESSION// /}.merged.fq

echo -e "fninishing merging paired end reads of $ACCESSION at $(date), results stored as $OUTPATH/$ACCESSION.merged.fq"

echo -e "start to align merged reads of $ACCESSION at $(date)"

$DIAMOND blastx  -d  $DIAMOND_DB  -q $OUTPATH/$ACCESSION.merged.fq -o $OUTPATH/$ACCESSION.merged.tab -f 6 -k 10 --threads 8 

echo -e "finishing alignment of $OUTPATH/$ACCESSION.merged.fq, start to genrate coverage at $(date)"

cat $OUTPATH/$ACCESSION.merged.tab | awk '{if($9 < $10){print($2"\t"$9-1"\t"$10)}else{print($2"\t"$10-1"\t"$9)}}' |  sort -k1,1 -k2,2n  > $OUTPATH/$ACCESSION.merged.bed 

cat $OUTPATH/$ACCESSION.merged.tab | awk '$3 >= 80 || $12 >= 120'  | awk '{if($9 < $10){print($2"\t"$9-1"\t"$10)}else{print($2"\t"$10-1"\t"$9)}}' |  sort -k1,1 -k2,2n > $OUTPATH/$ACCESSION.filtered.merged.bed


$BEDTOOLS  genomecov -i $OUTPATH/$ACCESSION.filtered.merged.bed -g $PROTLIST > $OUTPATH/$ACCESSION.filtered.merged.cov.txt 

rm $OUTPATH/$ACCESSION.merged.fq

echo -e "finishing alignment on $ACCESSION at $(date)"

done


echo -e "summarizing on SAMPLELIST at $(date)"

$RSCRIPT --vanilla $SUMUPR $ACCESSIONLIST $kvalue $phlist $OUTDIR1 $ABUNDANCEFILES $SIMULATIONDIR

echo -e "summarizing on SAMPLELIST completed at $(date)"



echo -e "job completed at $(date)"


INPATH="$OUTDIR1"
samplelist="$ACCESSIONLIST"

cat $samplelist | while read sample
do

file="$INPATH/$kvalue/$sample/$sample.merged.tab"
echo -e "kvalue\ttestset\tcluster\tgenome\tavrg_id_tp\tavrg_id_fp\tavrg_mm_tp\tavrg_mm_fp\tavrg_score_tp\tavrg_score_fp" > $INPATH/$kvalue/$sample/$sample.$kvalue.compstats.txt

cut -d "_" -f1,2 $file | sort | uniq  | grep -v "GUT" | while read gm
do
if [ $(grep -w "$gm" $phlist | wc -l) -gt 0 ]
then 

grep -w "$gm" $phlist |cut -d "|" -f4 | cut -d ":" -f1 | while read c
do

grep "${gm}_" $file | grep -w "$c"  > $INPATH/tmp1
grep "${gm}_" $file | grep -v "$c" > $INPATH/tmp2

#echo -e "$kvalue\t$sample\t$c\t$gm"

echo -e "$kvalue\t$sample\t$c\t$gm\t$(cat $INPATH/tmp1 | cut -f3 | awk '{SUM+=$1}END{print SUM / NR}')\t$(cat $INPATH/tmp2 | cut -f3 | awk '{SUM+=$1}END{print SUM / NR}')\t$(cat $INPATH/tmp1 | cut -f5 | awk '{SUM+=$1}END{print SUM / NR}')\t$(cat $INPATH/tmp2 | cut -f5 | awk '{SUM+=$1}END{print SUM / NR}')\t$(cat $INPATH/tmp1 | cut -f12 | awk '{SUM+=$1}END{print SUM / NR}')\t$(cat $INPATH/tmp2 | cut -f12 | awk '{SUM+=$1}END{print SUM / NR}')" 
done

fi
done >>  $INPATH/$kvalue/$sample/$sample.$kvalue.compstats.txt

echo -e "done with $file"
done
