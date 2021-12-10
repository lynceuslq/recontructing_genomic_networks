#!/bin/bash

kvalue="sensitive"
INPATH="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/abd_dist_model_testsets/onlyclusteredphage"
OUTDIR1="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/rm_overlap"
OUTDIR="$OUTDIR1/$kvalue"
REFORMATTER="/hwfssz5-tmp/ST_INFECTION/GlobalDatabase/user/fengqikai/software/bbmap/reformat.sh"
DIAMOND="/hwfssz5-tmp/ST_INFECTION/GlobalDatabase/user/fengqikai/software/Diamond/2.0.9/diamond"
DIAMOND_DB="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/rm_overlap/cl2_recoverprots.clustered.len.overlappingremoved.dmnd"
BEDTOOLS="/zfssz2/ST_MCHRI/COHORT/fengqikai/software/bedtools2/bedtools2/bin/bedtools"
PROTLIST="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/rm_overlap/cl2_recoverprots.clustered.len.overlappingremoved.length.txt"
ABUNDANCEFILES="$OUTDIR1/redoabd"
ACCESSIONLIST="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/rm_overlap/SAMPLELIST"
phlist="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/rm_overlap/gpdvcmorethan3.clustered.len.overlappingremoved.phage.list"
SIMULATIONDIR="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/abd_dist_model_testsets/abd_pack"
RSCRIPT="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/miniconda/bin/Rscript"
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

$DIAMOND blastx  -d  $DIAMOND_DB  -q $OUTPATH/$ACCESSION.merged.fq -o $OUTPATH/$ACCESSION.merged.tab -f 6 -k 10 --threads 32 

echo -e "finishing alignment of $OUTPATH/$ACCESSION.merged.fq, start to genrate coverage at $(date)"

#########cat $OUTPATH/$ACCESSION.merged.tab | awk '{if($9 < $10){print($2"\t"$9-1"\t"$10)}else{print($2"\t"$10-1"\t"$9)}}' |  sort -k1,1 -k2,2n  > $OUTPATH/$ACCESSION.merged.bed 

cat $OUTPATH/$ACCESSION.merged.tab | awk '$3 >= 85 && $12 >= 120'  | awk '{if($9 < $10){print($2"\t"$9-1"\t"$10)}else{print($2"\t"$10-1"\t"$9)}}' |  sort -k1,1 -k2,2n > $OUTPATH/$ACCESSION.filtered.merged.bed


$BEDTOOLS  genomecov -i $OUTPATH/$ACCESSION.filtered.merged.bed -g $PROTLIST > $OUTPATH/$ACCESSION.filtered.merged.cov.txt 

rm $OUTPATH/$ACCESSION.merged.fq

echo -e "finishing alignment on $ACCESSION at $(date)"

done


echo -e "summarizing on SAMPLELIST at $(date)"

$RSCRIPT --vanilla $SUMUPR $ACCESSIONLIST $kvalue $phlist $OUTDIR1 $ABUNDANCEFILES $SIMULATIONDIR

echo -e "summarizing on SAMPLELIST completed at $(date)"



echo -e "job completed at $(date)"
