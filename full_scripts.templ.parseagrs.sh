#!/bin/bash

SAMTOOLS="/ldfssz1/ST_META/share/User/zhujie/.conda/envs/bioenv/bin/samtools"
BOWTIE2="/ldfssz1/ST_META/share/User/zhujie/.conda/envs/bioenv/bin/bowtie2"
GPD_index="/ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/liqian6/fece_meta/GPD_index"
inputpath="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/abd_dist_model_testsets/clusteredwithhosts"
outputdir1="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/testongenomes/abd_dist_model_co_split"
samplelist="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/testongenomes/abd_dist_model_co_split/SAMPLIST"
genomelist="/ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/liqian6/fece_meta/genome.list"
ABUNDANCEFILES="$OUTDIR1/rabd_by_gm_cov10"
phlist="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/rm_overlap/gpdvcmorethan3.clustered.len.overlappingremoved.phage.list"
SIMULATIONDIR="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/abd_dist_model_testsets/abd_pack"

RSCRIPT="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/miniconda/bin/Rscript"
SUMUPR="/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/testongenomes/abd_dist_model_co_split/sum.gm.cov10.parseagrs.R"


while read sample

do

echo -e "mapping $sample to GPD genomes at $(date)" 
mkdir $outputdir1/${sample// /}
$BOWTIE2 -N 1 -p 32 -x $GPD_index -1 $inputpath/${sample// /}_R1.fastq -2 $inputpath/${sample// /}_R2.fastq -S $outputdir1/${sample// /}/${sample// /}.sam

outputdir="$outputdir1/${sample// /}"

echo -e "generating bam files for $sample"
$SAMTOOLS view -b -S -@ 32 $outputdir1/${sample// /}/${sample// /}.sam > $outputdir1/${sample// /}/${sample// /}.bam

echo -e "start sorting on  $sample" 
#sort a BAM file
cat $outputdir/${sample// /}.bam  | $SAMTOOLS sort -@ 4 -o $outputdir/${sample// /}_sorted.bam

#get comprehensive statistics
echo -e "start generating stats from  $sample" 
$SAMTOOLS stats -@ 8 $outputdir/${sample// /}_sorted.bam > $outputdir/${sample// /}_sorted.stats

#get coverage
echo -e "start generating indice for  $sample" 
$SAMTOOLS index -@ 8 $outputdir/${sample// /}_sorted.bam $outputdir/${sample// /}_sorted.bai

rm $outputdir1/${sample// /}/${sample// /}.sam 

#calculate average coverage
echo -e "start generating coverage on  $sample" 
/zfssz2/ST_MCHRI/COHORT/fengqikai/software/bedtools2/bedtools2/bin/bedtools bamtobed -i $outputdir/${sample// /}_sorted.bam > $outputdir/${sample// /}_sorted.bed


/zfssz2/ST_MCHRI/COHORT/fengqikai/software/bedtools2/bedtools2/bin/bedtools genomecov -i $outputdir/${sample// /}_sorted.bed  -g $genomelist  > $outputdir/${sample// /}_cov.txt

rm $outputdir/${sample// /}.bam

echo -e "$sample finished at $(date)"

rm $outputdir/${sample// /}_sorted.bam
done < $samplelist


mkdir $ABUNDANCEFILES 
echo -e "summarizing on $samplelist at $(date)"

$RSCRIPT --vanilla $SUMUPR $samplelist $phlist $outputdir1 $ABUNDANCEFILES $SIMULATIONDIR

echo -e "summarizing on $samplelist completed at $(date)"

echo -e "job completed"
date
