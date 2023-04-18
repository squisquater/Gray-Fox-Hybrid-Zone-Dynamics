#!/bin/bash -l
#SBATCH --job-name=bam2fastq
#SBATCH --array=1-x #number of samples
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --time 2-00:00:00
#SBATCH --mem=4GB
#SBATCH -p high
#SBATCH -o bam2fastq%A_%a.out
#SBATCH -e bam2fastq%A_%a.err

start=`date +%s`

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

REFERENCE=~/Reference_Genomes/canFam3_withY.fa #reference genome
SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" bamlist_depth.txt | cut -f1) #Sample ID
BAM=$(sed "${SLURM_ARRAY_TASK_ID}q;d" bamlist_depth.txt | cut -f2) #Path to subsetted bam file
MIN_DEPTH=$(sed "${SLURM_ARRAY_TASK_ID}q;d" bamlist_depth.txt | cut -f4) #minimum depth (0.5x sample mean)
MAX_DEPTH=$(sed "${SLURM_ARRAY_TASK_ID}q;d" bamlist_depth.txt | cut -f5) #maximum depth (2.5x sample mean)
OUT=/group/ctbrowngrp2/sophiepq/GrayFoxWGS/GrayFox1/angsd/PSMC/final_fqfiles #output directory

echo "converting .bam to .fastq for ${SAMPLE}"

module load htslib
module load samtools
module load bcftools

samtools mpileup -Q 30 -q 30 -u -v -f ${REFERENCE} ${BAM} |
bcftools call -c |  
vcfutils.pl vcf2fq -d ${MIN_DEPTH} -D ${MAX_DEPTH} -Q 30 | gzip > ${OUT}/${SAMPLE}.fq.gz

end=`date +%s`
runtime=$((end-start))
echo "finished converting ${SAMPLE}.bam to ${SAMPLE}.fastq after ${runtime} seconds"
