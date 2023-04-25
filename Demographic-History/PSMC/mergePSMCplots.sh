#!/usr/bin/env bash
#SBATCH --job-name=addlabel
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --array=1-50
#SBATCH --time 25:00:00
#SBATCH --mem=5GB
#SBATCH -o /group/ctbrowngrp2/sophiepq/GrayFoxWGS/GrayFox1/angsd/PSMC/PSMC_files_new/bootstrap/Plots/S19_3133/addlabel_%A_%a.out
#SBATCH -e /group/ctbrowngrp2/sophiepq/GrayFoxWGS/GrayFox1/angsd/PSMC/PSMC_files_new/bootstrap/Plots/S19_3133/addlabel_%A_%a.err
#SBATCH -p high
        
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
NUMBER=$(sed "${SLURM_ARRAY_TASK_ID}q;d" /group/ctbrowngrp2/sophiepq/GrayFoxWGS/GrayFox1/angsd/PSMC/PSMC_files_new/bootstrap/PSMC_S19_3133/Rounds.txt | cut -f1)

echo "58 ${NUMBER}" | awk '{ while(b<$1) {b=b+1; print $2}}' > new_${NUMBER}.txt ##this generates the run number and adds it to a new file 58 times because the output .txt files are 58 lines long.

paste S19_3133.${NUMBER}.txt new_${NUMBER}.txt > S19_3133.bootstrap_${NUMBER}.txt

### STEP THREE: CAT ALL THE FILES TOGETHER 

cat S19_3133.0.txt S19_3133.bootstrap* > S19_3133_allbootstrap_labeled.txt


# /group/ctbrowngrp2/sophiepq/GrayFoxWGS/GrayFox1/angsd/PSMC/PSMC_files_new/bootstrap/Plots/S19_3133/addlabel.sh
