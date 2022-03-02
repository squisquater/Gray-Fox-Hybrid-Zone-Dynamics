#!/bin/bash -l
#SBATCH --job-name=clean
#SBATCH --array=1-47
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --time 1-00:00:00
#SBATCH --mem=4GB
#SBATCH -p bmm
#SBATCH -A ctbrowngrp
#SBATCH -o /home/sophiepq/GrayFoxWGS/GrayFox1/slurmlogs/%A_%a.out
#SBATCH -e /home/sophiepq/GrayFoxWGS/GrayFox1/slurmlogs/%A_%a.err

start=`date +%s`

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

SID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ~/GrayFoxWGS/GrayFox1/align/mytrimmedfastqlist.txt | cut -f3)
SAMPLE=${SID}

echo "merging and cleaning up bams for ${SAMPLE}"

module load samtools

# merge 
samtools merge ~/GrayFoxWGS/GrayFox1/align/merged/${SAMPLE}.bam ~/GrayFoxWGS/GrayFox1/align/chunks/${SAMPLE}.*.bam

mid1=`date +%s`
runtime=$((mid1-start))
echo "finished merging chunks of ${SAMPLE} after ${runtime} seconds"

echo "changing read groups of ${SAMPLE}"
module load picardtools/2.7.1

java -jar ${PICARD}/picard.jar AddOrReplaceReadGroups \
      I=~/GrayFoxWGS/GrayFox1/align/${SAMPLE}.bam \
      O=~/GrayFoxWGS/GrayFox1/align/merged/RGreplaced/${SAMPLE}.bam \
      RGID=${SAMPLE} \
      RGLB=1 \
      RGPL=illumina \
      RGPU=NULL \
      RGSM=${SAMPLE}

mid2=`date +%s`
runtime=$((mid2-mid1))
echo "finished changing read groups of ${SAMPLE} after ${runtime} seconds"

# mark duplicates 
java -jar /home/hennelly/bin/picard.jar MarkDuplicates \
      I=~/GrayFoxWGS/GrayFox1/align/merged/RGreplaced/${SAMPLE}.bam \
      O=~/GrayFoxWGS/GrayFox1/align/merged/RGreplaced/dupmarked/${SAMPLE}.bam \
      M=~/GrayFoxWGS/GrayFox1/align/merged/RGreplaced/dupmarked/${SAMPLE}_metrics.txt

mid3=`date +%s`
runtime=$((mid3-mid2))
echo "finished removing duplicates from ${SAMPLE} after ${runtime} seconds"

# remove bad reads (duplicates, secondary alignments, mapQ<10)
samtools view -hb -F 256 -q 10 ~/GrayFoxWGS/GrayFox1/align/merged/RGreplaced/dupmarked/${SAMPLE}.bam | \
samtools view -hb -F 1024 > ~/GrayFoxWGS/GrayFox1/align/merged/RGreplaced/dupmarked/cleaned/${SAMPLE}_cleaned.bam

mid4=`date +%s`
runtime=$((mid4-mid3))
echo "finished removing bad reads from ${SAMPLE} after ${runtime} seconds"

# index
samtools index ~/GrayFoxWGS/GrayFox1/align/merged/RGreplaced/dupmarked/cleaned/${SAMPLE}_cleaned.bam

mid5=`date +%s`
runtime=$((mid5-mid4))
echo "finished indexing ${SAMPLE} after ${runtime} seconds"
end=`date +%s`
runtime=$((end-start))
echo "finished total cleaning for ${SAMPLE} after ${runtime} seconds"

################################################################################################
#################################Instructions for Running Script################################
################################################################################################
### --array = 1-n where n is the number of fastq files you want to trim.
### change file path to output (-o) and error (-e) files
### Note: You will need to include or exclude the merge command and update file paths accordingly 
### depending on whether you split your bam files up to speed up alignment.
### I don’t suggest mixing split and unsplit samples because the file paths won’t work for both 
### of them. 
### If you have some that are split up write a script for those and a different one for all the 
### unsplit .bam files
################################################################################################
