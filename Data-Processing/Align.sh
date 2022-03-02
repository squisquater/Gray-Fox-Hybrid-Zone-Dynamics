#!/bin/bash -l
#SBATCH --job-name=bwa
#SBATCH --array=1-47
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --time 5-00:00:00
#SBATCH --mem=7GB
#SBATCH -p bmm
#SBATCH -A ctbrowngrp
#SBATCH -o /home/sophiepq/GrayFoxWGS/GrayFox1/slurmlogs/%A_%a.out
#SBATCH -e /home/sophiepq/GrayFoxWGS/GrayFox1/slurmlogs/%A_%a.err

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

READ1=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ~/GrayFoxWGS/GrayFox1/align/mytrimmedfastqlist.txt | cut -f1)
READ2=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ~/GrayFoxWGS/GrayFox1/align/mytrimmedfastqlist.txt | cut -f2)
SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ~/GrayFoxWGS/GrayFox1/align/mytrimmedfastqlist.txt | cut -f3)
REF=~/Reference_Genomes/canFam3_withY.fa
INDIR=~/GrayFoxWGS/GrayFox1/trim/
OUTDIR=~/GrayFoxWGS/GrayFox1/align/

module load bwa
module load samtools

# Align R1 and R2 
# Pipe to a bam file that excludes bad mapping scores
# Pipe to a sorted bam for merging

cd $OUTDIR
bwa mem -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}" ${REF} ${INDIR}${READ1} ${INDIR}${READ2} | \
samtools view -q 10 -h -b - | \
samtools sort -o ${OUTDIR}${SAMPLE}.bam -


################################################################################################
#################################Instructions for Running Script################################
################################################################################################
### --array = 1-n where n is the number of trimmed fastq files you want to align.
### change file path to output (-o) and error (-e) files
### READ1 = List of names for trimmed R1 fastq files (e.g. ‘S13_2949_trim_R1.fastq.gz’)
### READ2 = List of names for trimmed R2 fastq files (e.g. ‘S13_2949_trim_R2.fastq.gz’)
### SAMPLE = List of the common SID name across both reads (e.g. ‘S13_2949_trim’)
### REF = Path to the reference genome
### INDIR = Path to directory with the trimmed fastq files
### OUTDIR = Path to the directory you want your aligned bam files to go into
################################################################################################
