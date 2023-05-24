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

## Local Ancestry Inference and Admixture Timing

## Geographic Cline Analaysis
#### Run cline models for: 
* All combined data [HZAR_ALLData.R]
* Western gray fox only [HZAR_West.R]
* Eastern gray fox only [HZAR_East.R]

## Demographic History 

## PSMC Analysis
**STEP0:** Subset autosomes [[subset_autosomes.sh]](https://github.com/squisquater/Gray-Fox-Hybrid-Zone-Dynamics/blob/main/Demographic-History/PSMC/subset_autosomes.sh) \
**STEP1:** Convert bam to fastq file [[bam2fastq.sh]](https://github.com/squisquater/Gray-Fox-Hybrid-Zone-Dynamics/blob/main/Demographic-History/PSMC/bam2fastq.sh) \
**STEP2:** Convert fastq to psmc fasta file [[fastq2psmcfa.sh]](https://github.com/squisquater/Gray-Fox-Hybrid-Zone-Dynamics/blob/main/Demographic-History/PSMC/fastq2psmcfa.sh) \
**STEP3:** Run PSMC for all files using Wang et al (2020) parameters [[runPSMC.sh]]() \
**STEP4:** Bootstrapping 
> Split fasta files for bootstrapping [[bootstrap_split.sh]](https://github.com/squisquater/Gray-Fox-Hybrid-Zone-Dynamics/blob/main/Demographic-History/PSMC/bootstrap_split.sh) \
> Run bootstrapping [[bootstrap.sh]](https://github.com/squisquater/Gray-Fox-Hybrid-Zone-Dynamics/blob/main/Demographic-History/PSMC/bootstrap.sh) 

**STEP5:** Merge all .psmc files together [[mergePSMC.sh]](https://github.com/squisquater/Gray-Fox-Hybrid-Zone-Dynamics/tree/main/Demographic-History/PSMC) \
**STEP6:** Plot PSMC trajectories [[plotPSMC.sh]]()

## Stairway Plot Analysis

## Tajima's D and Fu's F

## Heterozygosity

## Ecological Niche Modeling





