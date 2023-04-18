#!/bin/bash -l
#SBATCH --job-name=subset_autosomes
#SBATCH --array=1-x #number of samples
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --time 4:00:00
#SBATCH --mem=2GB
#SBATCH -p bmh
#SBATCH -A ctbrowngrp
#SBATCH -o /home/sophiepq/GrayFoxWGS/GrayFox1/slurmlogs/subset_autosomes_%a.out
#SBATCH -e /home/sophiepq/GrayFoxWGS/GrayFox1/slurmlogs/subset_autosomes_%a.err

SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" bamlist.txt | cut -f1) #sample ID
BAM=$(sed "${SLURM_ARRAY_TASK_ID}q;d" bamlist.txt | cut -f2) #path to bamefile
BEDIN= chr_keep.bed #bedfile with list of chromosomes to keep
BAMOUT= bam_subset_autosomes/ #directory to put subsetted bams in

module load bedtools

bedtools intersect -abam $BAM -b ${BEDIN} > ${BAMOUT}${SAMPLE}_subset.bam
