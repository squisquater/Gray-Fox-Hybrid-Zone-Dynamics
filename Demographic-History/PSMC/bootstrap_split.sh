#!/bin/bash -l
#SBATCH --job-name=split_psmcfa
#SBATCH --array=1-x #number of samples
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --time 12:00:00
#SBATCH --mem=4GB
#SBATCH -p high
#SBATCH -o split_psmcfa.out
#SBATCH -e split_psmcfa.err

start=`date +%s`

SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" sample_list.txt | cut -f1) #Sample ID
INFILE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" sample_list.txt | cut -f2) #Path to PSMC fasta file
OUTDIR=/group/ctbrowngrp2/sophiepq/GrayFoxWGS/GrayFox1/angsd/PSMC/PSMC_files_new

echo "splitting psmc.fa files for ${SAMPLE} for bootstrapping"

/group/ctbrowngrp2/hennelly/hennelly/bin/psmc/utils/splitfa ${INFILE} > ${OUTDIR}/${SAMPLE}_split.psmc.fa

end=`date +%s`
runtime=$((end-start))
echo "finished splitting psmc.fa for ${SAMPLE} bootstrapping after ${runtime} seconds"
