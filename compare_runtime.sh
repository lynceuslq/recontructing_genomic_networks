#!/bin/bash

kvalue="rep5"
INPATH="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/abd_dist_model_testsets/onlyclusteredphage"
OUTDIR1="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/testcomtime"
OUTDIR="$OUTDIR1/$kvalue"
REFORMATTER="/hwfssz5-tmp/ST_INFECTION/GlobalDatabase/user/fengqikai/software/bbmap/reformat.sh"
DIAMOND="/hwfssz5-tmp/ST_INFECTION/GlobalDatabase/user/fengqikai/software/Diamond/2.0.9/diamond"
DIAMOND_DB="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/rm_overlap/cl2_recoverprots.clustered.len.overlappingremoved.dmnd"
BEDTOOLS="/zfssz2/ST_MCHRI/COHORT/fengqikai/software/bedtools2/bedtools2/bin/bedtools"
PROTLIST="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/rm_overlap/cl2_recoverprots.clustered.len.overlappingremoved.length.txt"
ABUNDANCEFILES="$OUTDIR1/rep5abd"
ACCESSIONLIST="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/testcomtime/co.testtime.spl.ae"
phlist="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/rm_overlap/gpdvcmorethan3.clustered.len.overlappingremoved.phage.list"
SIMULATIONDIR="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/abd_dist_model_testsets/abd_pack"
RSCRIPT="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/miniconda/bin/Rscript"
SUMUPR="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/recoverprots_gpdvcmorethan3/sumup.argparse.protcovfilter.R"

mkdir $OUTDIR
mkdir $ABUNDANCEFILES

cat $ACCESSIONLIST | while read ACCESSION;
do
start=`date +%s`
echo -e "start to process fastq reads from $ACCESSION at $(date)" 
mkdir $OUTDIR/$ACCESSION

OUTPATH="$OUTDIR/$ACCESSION"

echo -e "start to merge paired end reads of $ACCESSION at $(date)"

$REFORMATTER in1=$INPATH/${ACCESSION// /}_R1.fastq  in2=$INPATH/${ACCESSION// /}_R2.fastq  out=$OUTPATH/${ACCESSION// /}.merged.fq

echo -e "fninishing merging paired end reads of $ACCESSION at $(date), results stored as $OUTPATH/$ACCESSION.merged.fq"

echo -e "start to align merged reads of $ACCESSION at $(date)"
start=`date +%s`

$DIAMOND blastx  -d  $DIAMOND_DB  -q $OUTPATH/$ACCESSION.merged.fq -o $OUTPATH/$ACCESSION.merged.tab -f 6 -k 10 --threads 1 

echo -e "finishing alignment of $OUTPATH/$ACCESSION.merged.fq, start to genrate coverage at $(date)"

cat $OUTPATH/$ACCESSION.merged.tab | awk '{if($9 < $10){print($2"\t"$9-1"\t"$10)}else{print($2"\t"$10-1"\t"$9)}}' |  sort -k1,1 -k2,2n  > $OUTPATH/$ACCESSION.merged.bed 

cat $OUTPATH/$ACCESSION.merged.tab | awk '$3 >= 80 || $12 >= 150'  | awk '{if($9 < $10){print($2"\t"$9-1"\t"$10)}else{print($2"\t"$10-1"\t"$9)}}' |  sort -k1,1 -k2,2n > $OUTPATH/$ACCESSION.filtered.merged.bed


$BEDTOOLS  genomecov -i $OUTPATH/$ACCESSION.filtered.merged.bed -g $PROTLIST > $OUTPATH/$ACCESSION.filtered.merged.cov.txt 

rm $OUTPATH/$ACCESSION.merged.fq

echo -e "finishing alignment on $ACCESSION at $(date)"

end=`date +%s`

runtime=$((end-start))
num=$(cat $INPATH/${ACCESSION// /}_R1.fastq | wc -l)
readnum=$(($num/2))

echo -e "$ACCESSION\t$readnum\tprot\t$start\t$end\t$runtime" >> $OUTDIR1/rep5.recordtime.txt

done
export PATH="/hwfssz5/ST_INFECTION/GlobalDatabase/user/fengqikai/software/.conda/envs/Trinity-2.11.0/bin/:$PATH"

GPD_index="/ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/liqian6/fece_meta/GPD_index"
inputpath="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/abd_dist_model_testsets/onlyclusteredphage"
outputdir1="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/testcomtime/rep5gm"
samplelist="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/testcomtime/co.testtime.spl.ae"
genomelist="/ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/liqian6/fece_meta/genome.list"
SAMTOOLS="/ldfssz1/ST_META/share/User/zhujie/.conda/envs/bioenv/bin/samtools"
BOWTIE2="/ldfssz1/ST_META/share/User/zhujie/.conda/envs/bioenv/bin/bowtie2"

while read sample

do

echo -e "mapping $sample to GPD genomes at $(date)" 
mkdir $outputdir1/${sample// /}

start=`date +%s`

$BOWTIE2 -N 1 -p 1 -x $GPD_index -1 $inputpath/${sample// /}_R1.fastq -2 $inputpath/${sample// /}_R2.fastq -S $outputdir1/${sample// /}/${sample// /}.sam

outputdir="$outputdir1/${sample// /}"

echo -e "generating bam files for $sample"
$SAMTOOLS view -b -S -@ 1 $outputdir1/${sample// /}/${sample// /}.sam > $outputdir1/${sample// /}/${sample// /}.bam

echo -e "start sorting on  $sample" 
#sort a BAM file
cat $outputdir/${sample// /}.bam  | $SAMTOOLS sort -@ 1 -o $outputdir/${sample// /}_sorted.bam

#get comprehensive statistics
echo -e "start generating stats from  $sample" 
$SAMTOOLS stats -@ 1 $outputdir/${sample// /}_sorted.bam > $outputdir/${sample// /}_sorted.stats

#get coverage
echo -e "start generating indice for  $sample" 
$SAMTOOLS index -@ 1 $outputdir/${sample// /}_sorted.bam $outputdir/${sample// /}_sorted.bai

rm $outputdir1/${sample// /}/${sample// /}.sam 

#calculate average coverage
echo -e "start generating coverage on  $sample" 
/zfssz2/ST_MCHRI/COHORT/fengqikai/software/bedtools2/bedtools2/bin/bedtools bamtobed -i $outputdir/${sample// /}_sorted.bam > $outputdir/${sample// /}_sorted.bed


/zfssz2/ST_MCHRI/COHORT/fengqikai/software/bedtools2/bedtools2/bin/bedtools genomecov -i $outputdir/${sample// /}_sorted.bed  -g $genomelist  > $outputdir/${sample// /}_cov.txt

rm $outputdir/${sample// /}.bam

echo -e "$sample finished at $(date)"

rm $outputdir/${sample// /}_sorted.bam

end=`date +%s`

runtime=$((end-start))
num=$(cat $inputpath/${sample// /}_R1.fastq | wc -l)
readnum=$(($num/2))

echo -e "$sample\t$readnum\tgenome\t$start\t$end\t$runtime" >> $outputdir1/rep5.recordtime.txt


done < $samplelist

echo -e "job completed"
date
