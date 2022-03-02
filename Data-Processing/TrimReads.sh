#!/bin/bash -l
#SBATCH --job-name=bbduk
#SBATCH --array=1-47
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -p bmm
#SBATCH -A ctbrowngrp
#SBATCH --mem 1GB 
#SBATCH -t 12:00:00
#SBATCH -o /home/sophiepq/GrayFoxWGS/GrayFox1/slurmlogs/01trim_%A_%a.out
#SBATCH -e /home/sophiepq/GrayFoxWGS/GrayFox1/slurmlogs/01trim_%A_%a.err

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

# Some tracking info
start=`date +%s`
echo $HOSTNAME
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

# Define variables
FASTQDIR=~/GrayFoxWGS/GrayFox1/tempfastq_untrimmed/
WORKDIR=~/GrayFoxWGS/GrayFox1/trim/
INNAME=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ~/GrayFoxWGS/GrayFox1/trim/myfastqtrimlist.txt | cut -f1)
SAMPNAME=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ~/GrayFoxWGS/GrayFox1/trim/myfastqtrimlist.txt | cut -f2)

module load bbmap
echo "Running BBduk on ${datapath}${sample}"
cd $WORKDIR

call="bbduk.sh -Xmx1g \
in1=${FASTQDIR}${INNAME}_R1_001.fastq.gz \
in2=${FASTQDIR}${INNAME}_R2_001.fastq.gz \
out1=${WORKDIR}${SAMPNAME}_trim_R1.fastq.gz \
out2=${WORKDIR}${SAMPNAME}_trim_R2.fastq.gz \
ref=~/GrayFoxWGS/GrayFox1/trim/adapters.fa \
ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=10"

echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo $runtime

################################################################################################
#################################Instructions for Running Script################################
################################################################################################
### --array = 1-n where n is the number of fastq files you want to trim.
### change file path to output (-o) and error (-e) files
### FASTQDIR = The directory where your fastq files are located
### WORKDIR = The location you want your trimmed reads to go
### INNAME = A list of file names (names include everything before the .fastq file designation)
### SAMPNAME = A list of desired output files names (just the SID)
################################################################################################
