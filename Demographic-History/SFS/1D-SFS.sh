#!/bin/bash -l
#SBATCH --job-name=1D-SFS
#SBATCH --nodes 1
#SBATCH --time 12:00:00
#SBATCH --mem 5GB
#SBATCH -p bmm
#SBATCH -A ctbrowngrp
#SBATCH -o /home/sophiepq/GrayFox/angsd/SFS/1D-SFS.out
#SBATCH -e /home/sophiepq/GrayFox/angsd/SFS/1D-SFS.err

module load angsd

echo "Generating 1-D Site Frequency Spectrums for all listed populations:"

cd /home/sophiepq/GrayFox/angsd/SFS

for POP in EastFullK3HC WestFullK3HC TexasFullK3HC CaliforniaFullK3HC
do
        echo $POP
        realSFS $POP.saf.idx > 1D.$POP.sfs
       
done
