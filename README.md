# Gray-Fox-Hybrid-Zone-Dynamics
Characterizing the dynamics of secondary contact between eastern and western gray fox lineages using whole genome sequencing

## **Data Processing**
* Trim reads [Trim.sh]
* Align to reference [Align.sh]
* Merge (if split) and clean reads [Merge_Clean.sh]
* Calculate read depth across samples [ReadDepth.sh]

## SNP Calling
* Genotype Likelihoods (angsd)

## Secondary Contact

### Local Ancestry Inference and Admixture Timing
#### Identify Ancestry Informative Markers
#### Remove Background LD
* Estimate distance based thinning parameter using reference population LD patterns [East_LD.sh] & [West_LD.sh]. Generate LD Decay Curves [LD_Curve.sh]
* Thin AIMs 
#### Convert physcal position (bp) to genetic position (cM)
#### Generate read counts (Maj/Min) for each admixed individual
#### Run AHMM

### Geographic Cline Analaysis

### Selection Across the Hybrid Zone


