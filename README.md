# Gray-Fox-Hybrid-Zone-Dynamics
Characterizing the dynamics of secondary contact between eastern and western gray fox lineages using whole genome and reduced-representation sequencing

## **Data Processing WGS**
* Trim reads [Trim.sh]
* Align to reference [Align.sh]
* Merge (if split) and clean reads [Merge_Clean.sh]
* Calculate read depth across samples [ReadDepth.sh]

## SNP Calling
* Genotype Likelihoods (angsd)

## **Data Processing GBS**
* See Kierepka et al. 2022 
* Details can be found in the github [repository](https://github.com/squisquater/Cryptic-Gray-Fox-Lineages-Secondary-Contact) associated with that manuscript.

## Demographic History 
### PSMC Analysis
**STEP0:** Subset autosomes [[subset_autosomes.sh]](https://github.com/squisquater/Gray-Fox-Hybrid-Zone-Dynamics/blob/main/Demographic-History/PSMC/subset_autosomes.sh) \
**STEP1:** Convert bam to fastq file [[bam2fastq.sh]](https://github.com/squisquater/Gray-Fox-Hybrid-Zone-Dynamics/blob/main/Demographic-History/PSMC/bam2fastq.sh) \
**STEP2:** Convert fastq to psmc fasta file [[fastq2psmcfa.sh]](https://github.com/squisquater/Gray-Fox-Hybrid-Zone-Dynamics/blob/main/Demographic-History/PSMC/fastq2psmcfa.sh) \
**STEP3:** Run PSMC for all files using Wang et al (2020) parameters [[runPSMC.sh]]() \
**STEP4:** Plot PSMC trajectories
 
### Stairway Plot Analysis
### Tajima's D and Fu's F

## Genetic Diversity
### Heterozygosity


## Secondary Contact

### Local Ancestry Inference and Admixture Timing
#### Identify Ancestry Informative Markers
#### Remove Background LD
* Estimate distance based thinning parameter using reference population LD patterns [East_LD.sh] & [West_LD.sh]. Generate LD Decay Curves [LD_Curve.sh]
* Thin AIMs 
#### Convert physcal position (bp) to genetic position (cM)
#### Generate read counts (Maj/Min) for each admixed individual
#### Run AHMM

see /group/ctbrowngrp2/sophiepq/GrayFoxWGS/GrayFox1/Ancestry_HMM/inputfiles_new/AIM0.6/IndividualAnalyses/LD_50kb for input and output files

### Geographic Cline Analaysis
#### Run cline models for: 
* All combined data [HZAR_ALLData.R]
* Western gray fox only [HZAR_West.R]
* Eastern gray fox only [HZAR_East.R]




