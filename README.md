# Gray-Fox-Hybrid-Zone-Dynamics
Characterizing the dynamics of secondary contact between eastern and western gray fox lineages using whole genome and reduced-representation sequencing

## **Data Processing WGS**
* Trim reads [Trim.sh]
* Align to reference [Align.sh]
* Merge (if split) and clean reads [Merge_Clean.sh]
* Calculate read depth across samples [ReadDepth.sh]

## **Data Processing GBS**
* See Kierepka et al. 2022 
* Details can be found in the github [repository](https://github.com/squisquater/Cryptic-Gray-Fox-Lineages-Secondary-Contact) associated with that manuscript.

### Local Ancestry Inference and Admixture Timing
#### Initial SNP Calling (GL)
```
#!/bin/bash -l
#SBATCH --job-name=angsd_allChrs_Q30_maf0.05_mind0.8_43GFs
#SBATCH --array=1-43
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --time 10:00:00
#SBATCH --mem=10GB
#SBATCH -p bmm
#SBATCH -A ctbrowngrp
#SBATCH -o /home/sophiepq/GrayFoxWGS/GrayFox1/slurmlogs/%A_%a.out
#SBATCH -e /home/sophiepq/GrayFoxWGS/GrayFox1/slurmlogs/%A_%a.err

OUTDIR=~/GrayFoxWGS/GrayFox1/angsd/vcf/
BAMLIST=~/GrayFoxWGS/GrayFox1/align/merged/RGreplaced/dupmarked/cleaned/GrayFox43.txt
REF=~/Reference_Genomes/canFam3_withY.fa
CHR=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ~/GrayFoxWGS/GrayFox1/angsd/chr_regions.txt | cut -f1)
MIN_DEPTH=100
MAX_DEPTH=430
MIN_IND=34

module load angsd

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

echo "Running angsd using dog-algined bams on region ${CHR}"

cd $OUTDIR

angsd -out ${OUTDIR}bychr_43GFs/angsd_${CHR}_Q30_maf0.05_mind0.8.vcf \
-r ${CHR} \
-ref ${REF} \
-bam ${BAMLIST} \
-remove_bads 1 \
-uniqueOnly 1 \
-minMapQ 30 -minQ 30 \
-doMajorMinor 2 \
-doCounts 1 -minMaf 0.05 -doMaf 8 \
-GL 1 -doGlf 2 \
-P 1 \
-setMinDepth ${MIN_DEPTH} \
-setMaxDepth ${MAX_DEPTH} \
-minInd ${MIN_IND}

echo "all done!"

# Analysis of genotype likelihoods (in ANGSD)

#I removed my failed individual (S16_1460) from my bam list to expedite the snp calling

# For calling variant sites in ANGSD, we do it chromosome by chromosome
# We use a non-HW dependent SNP calling algorithm based on counts of reads in ANGSD
# and filter out sites with unusually high or low coverage (> 2x or < .5x mean coverage)
# mean global coverage was ~225x (min = 110, max 450)
# then the samtools GL function to get genotype likelihoods for each individual,
# all filtering out low quality reads and bases
# Added in a minimum # individuals cutoff @ ~80%


# settings:
# -r specifies a genomic region to work on, e.g. chr1: specifies the entire chromosome
# -remove_bads removes reads with flags like duplicates
# -doMajorMinor 2: infer major and minor from allele counts
# -bam list of bams to include
# -GL 1: use samtools genotype likelihood method
# -doGlf 2: output beagle likelihood file
# -minMapQ 30 -minQ 20: filter out sites with low mapping quality or base/BAQ quality
# we use AlleleCounts method for MAF (-doMaf 8) with -doCounts
# which doesn't consider base quality and ignores all reads with non-target alleles
# but also doesn't rely on HW equlibrium assumptions to get ML allele freq from genotype freqs
# -minMaf x: and then do a cutoff to only include variant sites with >x minor allele freq.
# -minInd N: only keep sites with information (at least one read) from N individuals
# -P n means use n threads/nodes for each angsd task (here task=chromosome; then merges threads within-chrom)
# -setMaxDepth -setMinDepth filters out sites where total depth is below or exceeds some threshold
```
#### Merge all SNP-called chromosome files together
```
#!/bin/bash -l
#SBATCH --job-name=merge_allChr_GLfiles_43GFs
#SBATCH --time 1:00:00
#SBATCH --mem=1G
#SBATCH -p bmm
#SBATCH -A ctbrowngrp
#SBATCH -o /home/sophiepq/GrayFoxWGS/GrayFox1/slurmlogs/merge_allChr_GLfiles_43GFs.out
#SBATCH -e /home/sophiepq/GrayFoxWGS/GrayFox1/slurmlogs/merge_allChr_GLfiles_43GFs.err

DIR=~/GrayFoxWGS/GrayFox1/angsd/vcf/bychr_43GFs

## Start by making a temp file with just the header. You can pick any of the chr files to do this since you’re just pulling out the first line which should be the same for all of them.
#zcat ${DIR}/angsd_chr1:_Q30_maf0.05_mind0.8.vcf.beagle.gz | head -n 1 | gzip >> ${DIR}/allChr_header.gz

## Now pull all but the first line from every chromosome file and gzip it to a new file
for i in ${DIR}/angsd_chr*:_Q30_maf0.05_mind0.8.vcf.beagle.gz ; do
    zcat $i | tail -n +2 | gzip >> ${DIR}/allChr_GLfiles.beagle.gz
done

# Then cat the header file and the merged chromosome file together.
#zcat ${DIR}/allChr_header.gz ${DIR}/allChr_GLfiles.beagle.gz | gzip >> ${DIR}/allchr_GLfiles_header.beagle.gz
```

#### Calculate allele frequencies by calling genotypes in parental pops using the above determined SNP set.
**1. Make a SNP list to call genotypes in my parental reference populations.**
```
### Can run this interactively. It goes pretty fast
### Extract the first few column out of the GL file which has the SNP data
zcat ~/GrayFoxWGS/GrayFox1/angsd/vcf/bychr_43GFs/allchr_GLfiles_header.beagle.gz | cut -f 1-3 | gzip >> All_SNP_POS.txt

### Remove the header and split it up into two columns
zcat All_SNP_POS.txt | tail -n +2 | tr _ $'\t' | gzip >> 19.7m_Snp_Set_ForGLs.txt

##Extract the Major and Minor allele columns from the GL file and remove the header
zcat ~/GrayFoxWGS/GrayFox1/angsd/vcf/bychr_43GFs/allchr_GLfiles_header.beagle.gz | cut -f2,3 | gzip >> All_SNP_MAJMIN.txt
zcat ~/GrayFoxWGS/GrayFox1/angsd/vcf/bychr_43GFs/All_SNP_MAJMIN.txt | tail -n +2 | gzip >> All_SNP_MAJMIN_nohead.txt

Merge the two files back together
paste <(zcat 19.7m_Snp_Set_ForGLs.txt) <(zcat All_SNP_MAJMIN_nohead.txt) | gzip > 19.7m_Snp_Set_MajMin_ForGLs.txt

### Index/sort the file
module load angsd
angsd sites index 19.7m_Snp_Set_MajMin_ForGLs.txt
```
**2. Call genotypes in parental pops** 
```
#!/bin/bash -l
#SBATCH --job-name=call_eastRef_genotypes_14GFs
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --time 2-00:00:00
#SBATCH --mem=30GB
#SBATCH -p bmm
#SBATCH -A ctbrowngrp
#SBATCH -o /home/sophiepq/GrayFoxWGS/GrayFox1/slurmlogs/call_eastRef_genotypes_14GFs.out
#SBATCH -e /home/sophiepq/GrayFoxWGS/GrayFox1/slurmlogs/call_eastRef_genotypes_14GFs.err

POP=East_ref_14GFs
BAM_LIST=~/GrayFoxWGS/GrayFox1/align/merged/RGreplaced/dupmarked/cleaned/${POP}.txt
REF=~/Reference_Genomes/canFam3_withY.fa
SNP_FILE=~/GrayFoxWGS/GrayFox1/angsd/vcf/bychr_43GFs/19.7m_Snp_Set_MajMin_ForGLs.txt
OUTDIR=~/GrayFoxWGS/GrayFox1/angsd/vcf/

module load angsd

echo "Calculating allele frequencies for POP: ${POP}"

echo "finding site allele frequencies"
angsd -out ${OUTDIR}bychr_43GFs/angsd_${POP}_genotypes_19.7mSNPs \
-bam ${BAM_LIST} \
-ref ${REF} \
-underFlowProtect 1 \
-remove_bads 1 \
-doGeno 2 \
-doPost 1 \
-geno_minDepth 1 \
-geno_minMM 0.6 \
-minMapQ 30 \
-doMajorMinor 3 \
-doMaf 1 \
-GL 1 \
-sites ${SNP_FILE} \
-doCounts 1 \
-P 1

echo "done getting allele frequencies!"

# options
# basic quality filtering for reads
# -rf limits angsd to walking through indexed regions in the regions file, not the whole genome
# -doGeno 2 writes genotypes as 0, 1, 2 copies of the minor allele (or -1 for no genotype called)
# -doPost 1 estimates genotype posteriors using allele freq and HW as prior
# -geno_minDepth=6 requires a minimum depth of 6 reads to call a genotype.
# -geno_minMM=0.6 requires at minimum 60% of reads to match the major or minor allele.
# -GL 1 uses samtools GL method
# -doMajorMinor 3: takes major & minor allele from sites file
# ANGSD calculates freq of minor allele (ignoring all other non maj-min alleles seen)
# -doMaf 1 uses known major and minor alleles + gets a ML estimate of the allele freq based on EM and the genotype likelihoods (note: assumes HW). http://www.popgen.dk/angsd/index.php/Allele_Frequencies
# -minMapQ 30: filter out sites with low mapping quality. No filter for base/BAQ quality (error rates are incorporated into the model)
# (I pre-computed BAQ scores and replaced quality with minimum of BAQ/base quality,
# so this is equivalend to -baq 2 option here)
``` 

```
#!/bin/bash -l
#SBATCH --job-name=call_westRef_genotypes_14GFs
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --time 2-00:00:00
#SBATCH --mem=30GB
#SBATCH -p bmm
#SBATCH -A ctbrowngrp
#SBATCH -o /home/sophiepq/GrayFoxWGS/GrayFox1/slurmlogs/call_westRef_genotypes_14GFs.out
#SBATCH -e /home/sophiepq/GrayFoxWGS/GrayFox1/slurmlogs/call_westRef_genotypes_14GFs.err

POP=West_ref_14GFs
BAM_LIST=~/GrayFoxWGS/GrayFox1/align/merged/RGreplaced/dupmarked/cleaned/${POP}.txt
REF=~/Reference_Genomes/canFam3_withY.fa
SNP_FILE=~/GrayFoxWGS/GrayFox1/angsd/vcf/bychr_43GFs/19.7m_Snp_Set_MajMin_ForGLs.txt
OUTDIR=~/GrayFoxWGS/GrayFox1/angsd/vcf/

module load angsd

echo "Calculating allele frequencies for POP: ${POP}"

echo "finding site allele frequencies"
angsd -out ${OUTDIR}bychr_43GFs/angsd_${POP}_genotypes_19.7mSNPs \
-bam ${BAM_LIST} \
-ref ${REF} \
-underFlowProtect 1 \
-remove_bads 1 \
-doGeno 2 \
-doPost 1 \
-geno_minDepth 1 \
-geno_minMM 0.6 \
-minMapQ 30 \
-doMajorMinor 3 \
-doMaf 1 \
-GL 1 \
-sites ${SNP_FILE} \
-doCounts 1 \
-P 1

echo "done getting allele frequencies!"

# options
# basic quality filtering for reads
# -rf limits angsd to walking through indexed regions in the regions file, not the whole genome
# -doGeno 2 writes genotypes as 0, 1, 2 copies of the minor allele (or -1 for no genotype called)
# -doPost 1 estimates genotype posteriors using allele freq and HW as prior
# -geno_minDepth=6 requires a minimum depth of 6 reads to call a genotype.
# -geno_minMM=0.6 requires at minimum 60% of reads to match the major or minor allele.
# -GL 1 uses samtools GL method
# -doMajorMinor 3: takes major & minor allele from sites file
# ANGSD calculates freq of minor allele (ignoring all other non maj-min alleles seen)
# -doMaf 1 uses known major and minor alleles + gets a ML estimate of the allele freq based on EM and the genotype likelihoods (note: assumes HW). http://www.popgen.dk/angsd/index.php/Allele_Frequencies
# -minMapQ 30: filter out sites with low mapping quality. No filter for base/BAQ quality (error rates are incorporated into the model)
# (I pre-computed BAQ scores and replaced quality with minimum of BAQ/base quality,
# so this is equivalend to -baq 2 option here)
```
***3. Generate Ancestry Informative Markers***
*In R*
```
library('dplyr')
library('stringr')

setwd("C:/Users/Sophie/Desktop/GrayFox/WGS/angsd")

#read in files
EastRef <- read.table("angsd_East_ref_14GFs_genotypes_19.7mSNPs_genoD1.mafs.txt", header = TRUE)
WestRef <- read.table("angsd_West_ref_14GFs_genotypes_19.7mSNPs_genoD1.mafs.txt", header = TRUE)

#Make a new column in each that has the chr & pos columns concatenated
EastRef$chromo_pos <- str_c(EastRef$chromo,'_',EastRef$position)
WestRef$chromo_pos <- str_c(WestRef$chromo,'_',WestRef$position)

#merge two datasets - match by chr_pos column
EastWest_merged <- merge(x = EastRef, y = WestRef, by = "chromo_pos", all = TRUE)
#sapply(east_west_merged, class)

##subset data where the minimum number of individuals is set to 12 for both the eastern and western pops
Eastmin12 <- filter(EastWest_merged, nInd.x > 11)
Eastmin12_Westmin12 <- filter(Eastmin6, nInd.y >11)

##create a new column that calculates the AF difference between the two pops
Eastmin10_Westmin12$AF_Diff <- abs(Eastmin12_Westmin12$knownEM.x - Eastmin12_Westmin12$knownEM.y)

##Filter SNPS by AF_Diff at different thresholds
#EastWest_AFDiff0.2_ref12 <- filter(Eastmin12_Westmin12, AF_Diff > 0.2 )
#EastWest_AFDiff0.3_ref12 <- filter(Eastmin12_Westmin12, AF_Diff > 0.3 )
#EastWest_AFDiff0.4_ref12 <- filter(Eastmin12_Westmin12, AF_Diff > 0.4 )
#EastWest_AFDiff0.5_ref12 <- filter(Eastmin12_Westmin12, AF_Diff > 0.5 )
EastWest_AFDiff0.6_ref12 <- filter(Eastmin12_Westmin12, AF_Diff > 0.6 )
#EastWest_AFDiff0.7_ref12 <- filter(Eastmin12_Westmin12, AF_Diff > 0.7 )
#EastWest_AFDiff0.8_ref12 <- filter(Eastmin12_Westmin12, AF_Diff > 0.8 )
#EastWest_AFDiff0.9_ref12 <- filter(Eastmin12_Westmin12, AF_Diff > 0.9 )

#Write Files
#write.table(EastWest_AFDiff0.2, "EastWest_AFDiff0.2_genodepth1.txt")
#write.table(EastWest_AFDiff0.3, "EastWest_AFDiff0.3_genodepth1.txt")
#write.table(EastWest_AFDiff0.4, "EastWest_AFDiff0.4_genodepth1.txt")
#write.table(EastWest_AFDiff0.5, "EastWest_AFDiff0.5_genodepth1.txt")
write.table(EastWest_AFDiff0.6, "EastWest_AFDiff0.6_genodepth1.txt")
#write.table(EastWest_AFDiff0.7, "EastWest_AFDiff0.7_genodepth1.txt")
#write.table(EastWest_AFDiff0.8, "EastWest_AFDiff0.8_genodepth1.txt")
#write.table(EastWest_AFDiff0.9, "EastWest_AFDiff0.9_genodepth1.txt")
```
**4. Remove Background LD**
**Generate LD Decay Curves in Reference Pops**
**Input Files**
  * allchr_eastRef.beagle.gz
  * allchr_westRef.beagle.gz
 ```
 ~/bin/ngsTools/ngsLD/ngsLD --n_threads 10 --verbose 1 --geno ~/GrayFoxWGS/GrayFox1/angsd/vcf/bychr_43GFs/allchr_eastRef.beagle.gz --probs yes --max_kb_dist 2000 --rnd_sample 0.01 --pos ~/GrayFoxWGS/GrayFox1/angsd/vcf/bychr_43GFs/19.7m_Snp_Set_ForGLs.txt --n_ind 14 --n_sites 19760937 --out ~/GrayFoxWGS/GrayFox1/angsd/ngsLD/outputs/allchr_EastRefLD_thin0.01.id
 
 ~/bin/ngsTools/ngsLD/ngsLD --n_threads 10 --verbose 1 --geno ~/GrayFoxWGS/GrayFox1/angsd/vcf/bychr_43GFs/allchr_eastRef.beagle.gz --probs yes --max_kb_dist 2000 --rnd_sample 0.01 --pos ~/GrayFoxWGS/GrayFox1/angsd/vcf/bychr_43GFs/19.7m_Snp_Set_ForGLs.txt --n_ind 14 --n_sites 19760937 --out ~/GrayFoxWGS/GrayFox1/angsd/ngsLD/outputs/allchr_EastRefLD_thin0.01.id
```
**Plot LD Decay Curves**
```
Rscript --vanilla --slave ~/bin/ngsTools/ngsLD/scripts/fit_LDdecay.R --ld_files ~/GrayFoxWGS/GrayFox1/angsd/ngsLD/outputs/allchr_EastWest_0.001_list.txt --plot_x_lim 1000 --fit_level 300 --plot_scale 3 --out ~/GrayFoxWGS/GrayFox1/angsd/ngsLD/outputs/allchr_EastWestLD.pdf
```
where *allchr_EastWest_0.001_list.txt* is a .txt file that has a list of the two LD datasets (.id files with filepaths) you generated above

 **2. Thin Sites Based on 50kb Threshold**
 ```
 library(tidyverse)

setwd("C:/Users/Sophie/Desktop/GrayFox/WGS/angsd/")

AIM <- 0.6 ###Change Me!

bigfile <- read.table(paste0("EastWest_AFDiff", AIM, "_genodepth1_12ind.txt", header=TRUE)
chr <- c(paste0("chr", 1:38), "chrX")
windows <- c(25000, 50000, 75000, 100000)

for (i in 1:length(chr)){
      CHR <- chr[i]
      littlefile <- bigfile %>%
        filter(chromo.x == CHR) %>%
        arrange(position.x)
      
      # blank dataframe for True or False
      window_df <- as.data.frame(matrix(ncol = length(windows), nrow=nrow(littlefile)))
      names(window_df) <- paste0(windows/1000,"KB")

        for (j in 1:length(windows)){
          
          spacing <- windows[j]
          keep_snp <- rep(F, nrow(littlefile)) 
          last_rpos <- range(littlefile$position.x)[2] * -1
          for (i in 1:nrow(littlefile)) { 
            this_rpos <- littlefile[i, "position.x"] 
            if (this_rpos - last_rpos >= spacing){
              keep_snp[i] <- T 
              last_rpos <- this_rpos 
          } #else do nothing  
            window_df[,j] <- keep_snp
        }
        }
      newlittlefile <- bind_cols(littlefile, window_df)
      write.table(newlittlefile, paste0("AIM",AIM,"_",CHR,"_summary.txt"), sep="\t", quote = FALSE, col.names = TRUE)
}
```


#### Convert physcal position (bp) to genetic position (cM)
#### Generate read counts (Maj/Min) for each admixed individual
#### Run AHMM

## Demographic History 
### PSMC Analysis
**STEP0:** Subset autosomes [[subset_autosomes.sh]](https://github.com/squisquater/Gray-Fox-Hybrid-Zone-Dynamics/blob/main/Demographic-History/PSMC/subset_autosomes.sh) \
**STEP1:** Convert bam to fastq file [[bam2fastq.sh]](https://github.com/squisquater/Gray-Fox-Hybrid-Zone-Dynamics/blob/main/Demographic-History/PSMC/bam2fastq.sh) \
**STEP2:** Convert fastq to psmc fasta file [[fastq2psmcfa.sh]](https://github.com/squisquater/Gray-Fox-Hybrid-Zone-Dynamics/blob/main/Demographic-History/PSMC/fastq2psmcfa.sh) \
**STEP3:** Run PSMC for all files using Wang et al (2020) parameters [[runPSMC.sh]]() \
**STEP4:** Bootstrapping 
> Split fasta files for bootstrapping [[bootstrap_split.sh]](https://github.com/squisquater/Gray-Fox-Hybrid-Zone-Dynamics/blob/main/Demographic-History/PSMC/bootstrap_split.sh) \
> Run bootstrapping [[bootstrap.sh]](https://github.com/squisquater/Gray-Fox-Hybrid-Zone-Dynamics/blob/main/Demographic-History/PSMC/bootstrap.sh) 

**STEP5:** Merge all .psmc files together [[mergePSMC.sh]](https://github.com/squisquater/Gray-Fox-Hybrid-Zone-Dynamics/tree/main/Demographic-History/PSMC) \
**STEP6:** Plot PSMC trajectories [[plotPSMC.sh]]()

### Stairway Plot Analysis
### Tajima's D and Fu's F

## Genetic Diversity
### Heterozygosity

## Ecological Niche Modeling

## Secondary Contact

### Geographic Cline Analaysis
#### Run cline models for: 
* All combined data [HZAR_ALLData.R]
* Western gray fox only [HZAR_West.R]
* Eastern gray fox only [HZAR_East.R]




