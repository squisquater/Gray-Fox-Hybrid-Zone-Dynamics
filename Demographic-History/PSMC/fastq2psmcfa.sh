#!/bin/bash -l
#SBATCH --job-name=fq2psmcfa
#SBATCH --array=1-x #number of samples
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --time 6:00:00
#SBATCH --mem=4GB
#SBATCH -p high
#SBATCH -o fq2psmcfa%A_%a.out
#SBATCH -e fq2psmcfa%A_%a.err

start=`date +%s`

echo "converting .fastq files to psmc.fa files"

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" fqlist.txt | cut -f1) #Sample ID
FQ=$(sed "${SLURM_ARRAY_TASK_ID}q;d" fqlist.txt | cut -f2) #Path to fastq file
OUT=~/GrayFoxWGS/GrayFox1/angsd/PSMC/PSMC_files_new

/home/hennelly/bin/psmc/utils/fq2psmcfa ${FQ} > ${OUT}/${SAMPLE}.psmc.fa

end=`date +%s`
runtime=$((end-start))
echo "finished converting .fastq files after ${runtime} seconds"
