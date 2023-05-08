#!/bin/bash -l
#SBATCH --job-name=EastPureK3HC_doSaf
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --time 10:00:00
#SBATCH --mem=4GB
#SBATCH -p high
#SBATCH -o /home/sophiepq/GrayFox/angsd/SFS/EastPureK3HC_doSaf.out
#SBATCH -e /home/sophiepq/GrayFox/angsd/SFS/EastPureK3HC_doSaf.err

module load angsd

echo "Generating site allele frequency file using GLs called using samtools (-GL 1)"

POP=EastPureK3HC
REF=/home/sophiepq/Reference_Genomes/canFam3_withY.fa
ANC=/home/sophiepq/GrayFox/angsd/SFS/AncestralState/CanVulpAnc.fa
WORKDIR=/home/sophiepq/GrayFox/angsd/SFS/
OUT=${POP}
BAMLIST=/home/sophiepq/GrayFox/angsd/SFS/${POP}_bamlist.txt
MIN_DEPTH=483 # [0.5 * Pop Global Mean]
MAX_DEPTH=2415 # [2.5 * Pop Global Mean]
MIN_IND=10

# set filters
FILTERS="-minMapQ 20 -minQ 20 -baq 1 -uniqueOnly 1 -remove_bads 1 -trim 0 -setMinDepth ${MIN_DEPTH} -setMaxDepth ${MAX_DEPTH} -minInd ${MIN_IND} -doCounts 1"
OPT=" -dosaf 1 -gl 1"

cd $WORKDIR

# get site frequency likelihoods
angsd -P 2 -out $OUT -b $BAMLIST -anc $ANC -ref $REF $FILTERS $OPT
