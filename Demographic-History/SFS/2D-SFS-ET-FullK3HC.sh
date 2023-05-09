#!/bin/bash -l
#SBATCH --job-name=2D-SFS-ET-FullK3HC
#SBATCH --nodes 1
#SBATCH --time 12:00:00
#SBATCH --mem 5GB
#SBATCH -p bmm
#SBATCH -A ctbrowngrp
#SBATCH -o /home/sophiepq/GrayFox/angsd/SFS/2D-SFS-ET-FullK3HC.out
#SBATCH -e /home/sophiepq/GrayFox/angsd/SFS/2D-SFS-ET-FullK3HC.err

module load angsd

echo "Generating 2-D Site Frequency Spectrum for all listed populations:"

cd /home/sophiepq/GrayFox/angsd/SFS

POP0=EastFullK3HC #n=69 > 138 gene copies
POP1=TexasFullK3HC #n=87 > 174 gene copies

realSFS ${POP0}.saf.idx ${POP1}.saf.idx > 2D-SFS-ET-FullK3HC.sfs
