#!/bin/bash -l
#SBATCH --job-name=runPSMC
#SBATCH --array=1-x
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --time 12:00:00
#SBATCH --mem=2GB
#SBATCH -p high
#SBATCH -o runPSMC_%A_%a.out
#SBATCH -e runPSMC_%A_%a.err

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
STARTTIME=$(date +"%s")

echo "running PSMC analysis on $SAMPLE using the atomic intervals from Wang et al 2020"

SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" sample.txt | cut -f1) #sample ID
IN=~/GrayFoxWGS/GrayFox1/angsd/PSMC/PSMC_files_new
OUT=~/GrayFoxWGS/GrayFox1/angsd/PSMC/PSMC_files_new/PSMC_t15_pwang_r5

/home/hennelly/bin/psmc/psmc -t15 -r5 -p "4+25*2+4+6" -o  ${OUT}/${SAMPLE}_t15_pwang_r5.psmc ${IN}/${SAMPLE}.psmc.fa

ENDTIME=$(date +%s)
TIMESPEND=$(($ENDTIME - $STARTTIME))
((sec=TIMESPEND%60,TIMESPEND/=60, min=TIMESPEND%60, hrs=TIMESPEND/60))
timestamp=$(printf "%d:%02d:%02d" $hrs $min $sec)
echo "Took $timestamp hours:minutes:seconds to complete..."
