This folder contains information for generating a site frequency spectrum which is used in multiple downstream analyses including:
* Stairwayplots
* Tajima's D
* Heterozygosity

#### Target Populations Include:

#### Several datasets

--- fastsimcoal2 Dataset ---
##### Full K=3 Dataset High Coverage (n=182)
*I think this one will be the best representation for the fastsimcoal modeling*
* East All (n=69) **Bamlist:** ~/GrayFox/angsd/SFS/EastFullK3HC_bamlist.txt **DONE**
* West All (n=113) **Bamlist:** ~/GrayFox/angsd/SFS/WestFullK3HC_bamlist.txt **DONE**
  * Texas All (n=87) **Bamlist:** ~/GrayFox/angsd/SFS/TexasFullK3HC_bamlist.txt **DONE**
  * California All (n=26) **Bamlist:** ~/GrayFox/angsd/SFS/CaliforniaFullK3HC_bamlist.txt **DONE**

--- StairwayplotDatasets ---
Does the high coverage data look similar to the all coverage data for Tajima's D and Stairwayplots?
##### Pure Only Dataset All Coverage (n=195)
* East All (n=81) **Bamlist:** ~/GrayFox/angsd/SFS/EastPureAll_bamlist.txt **DONE**
* West All (n=114) **Bamlist:** ~/GrayFox/angsd/SFS/WestPureAll_bamlist.txt **DONE**
  * Texas All (n=42) **Bamlist:** ~/GrayFox/angsd/SFS/TexasK3PureAll_bamlist.txt **DONE**
  * California All (n=19) **Bamlist:** ~/GrayFox/angsd/SFS/AllCaliforniaPureAll_bamlist.txt **DONE**
  * N. California All (n=11) **Bamlist:** ~/GrayFox/angsd/SFS/N.CAPureAll_bamlist.txt **DONE**

##### Pure Only High Coverage Dataset
* East High Coverage (n=61) **Bamlist:** ~/GrayFox/angsd/SFS/East_bamlist.txt **DONE**
* West High Coverage (n=86) **Bamlist:** ~/GrayFox/angsd/SFS/West_bamlist.txt **DONE**
* Texas High Coverage (n=27) **Bamlist:** ~/GrayFox/angsd/SFS/Texas_bamlist.txt **DONE**
* California High Coverage (n=13) **Bamlist:** ~/GrayFox/angsd/SFS/CA_merged_bamlist.txt
* N. California High Coverage (n=5) **Bamlist:** ~/GrayFox/angsd/SFS/N.CA_bamlist.txt **DONE**


##### High Coverage K2 Conservative Dataset
* East K2Conservative High Coverage (n=60) **Bamlist:** 
* West K2Conservative High Coverage (n=40) **Bamlist:** 
* TexasK2Conservative High Coverage (n=0) *All Texas Individuals have some eastern introgression at K=2
* California High Coverage (n=13) **Bamlist:** ~/GrayFox/angsd/SFS/CA_merged_bamlist.txt
* N. California High Coverage (n=5) **Bamlist:** ~/GrayFox/angsd/SFS/N.CA_bamlist.txt

### STEP00: Calculate Depth for all samples --> [getdepths.sh]()
/home/sophiepq/GrayFox/angsd/SFS/HeterozygosityCalcs/getdepths.sh

### STEP0: Generate an ancestral fasta file to generate an unfolded SFS.
Include: [golden jackal + coyote + gray wolf] & [red fox + kit fox + arctic fox] \
Script: [CanVulpAncestralFasta.sh]()

### STEP1: Find a ‘global estimate’ of the SFS
#### STEP1a: Estimate the site allele frequency likelihood
*Note: If you don't have the ancestral state, you can instead estimate the folded SFS. This is done by supplying the -anc with the reference genome and applying -fold 1 to realSFS. (aka do not supply a -ref)*

[/home/sophiepq/GrayFox/angsd/SFS/EastFullK3HC_doSaf.sh]() \
[/home/sophiepq/GrayFox/angsd/SFS/WestFullK3HC_doSaf.sh]() \
[/home/sophiepq/GrayFox/angsd/SFS/TexasFullK3HC_doSaf.sh]() \
[/home/sophiepq/GrayFox/angsd/SFS/CaliforniaFullK3HC_doSaf.sh]() 

#### STEP1b: Generate a 1-D site frequency spectrum --> [1D-SFS.sh](https://github.com/squisquater/Gray-Fox-Hybrid-Zone-Dynamics/blob/main/Demographic-History/SFS/1D-SFS.sh)

