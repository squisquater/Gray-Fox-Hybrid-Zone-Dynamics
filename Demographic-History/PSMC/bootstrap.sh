#!/usr/bin/env bash
#SBATCH --job-name=S19_3133_bootstrap
#SBATCH --array=1-100
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --time 6:00:00
#SBATCH --mem=4GB
#SBATCH -o S19_3133_bootstrap_%A_%a.out
#SBATCH -e S19_3133_bootstrap_%A_%a.err
#SBATCH -p med

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
STARTTIME=$(date +"%s")

SAMPLE=S19_3133
ROUND=$(sed "${SLURM_ARRAY_TASK_ID}q;d" /group/ctbrowngrp2/sophiepq/GrayFoxWGS/GrayFox1/angsd/PSMC/PSMC_files_new/bootstrap/PSMC_S19_3133/Rounds.txt | cut -f1) #just a file with a list from 1-100. You could definitely set this script up differently

echo "running run ${ROUND} of the PSMC analysis bootstrapping on ${SAMPLE} using the atomic intervals from Wang et al 2020"

cd /group/ctbrowngrp2/sophiepq/GrayFoxWGS/GrayFox1/angsd/PSMC/PSMC_files_new/bootstrap/PSMC_${SAMPLE}/

seq 1 | xargs -i echo /group/ctbrowngrp2/hennelly/hennelly/bin/psmc/psmc -t15 -r5 -b -p "4+25*2+4+6" \
	   -o run-${ROUND}.psmc /group/ctbrowngrp2/sophiepq/GrayFoxWGS/GrayFox1/angsd/PSMC/PSMC_files_new/bootstrap/PSMC_${SAMPLE}/${SAMPLE}_autosomes_split.psmc.fa | sh

ENDTIME=$(date +%s)
TIMESPEND=$(($ENDTIME - $STARTTIME))
((sec=TIMESPEND%60,TIMESPEND/=60, min=TIMESPEND%60, hrs=TIMESPEND/60))
timestamp=$(printf "%d:%02d:%02d" $hrs $min $sec)
echo "Took $timestamp hours:minutes:seconds to complete..."
