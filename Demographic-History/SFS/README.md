This folder contains information for generating a site frequency spectrum which is used in multiple downstream analyses including:
* Stairwayplots
* Tajima's D
* Heterozygosity

#### Target Populations Include:
##### Full K=3 Dataset High Coverage
* East All (n=69) **Bamlist:** ~/GrayFox/angsd/SFS/EastFullK3HC_bamlist.txt **DONE**
* Texas All (n=87) **Bamlist:** ~/GrayFox/angsd/SFS/TexasFullK3HC_bamlist.txt **DONE**

##### Pure Only High Coverage Dataset
* East High Coverage (n=61) **Bamlist:** ~/GrayFox/angsd/SFS/East_bamlist.txt **DONE**
* Texas High Coverage (n=27) **Bamlist:** ~/GrayFox/angsd/SFS/Texas_bamlist.txt **DONE**

### STEP00: Calculate Depth for all samples --> [getdepths.sh]()
/home/sophiepq/GrayFox/angsd/SFS/HeterozygosityCalcs/getdepths.sh

### STEP0: Generate an ancestral fasta file to generate an unfolded SFS.
Include: [golden jackal + coyote + gray wolf] & [red fox + kit fox + arctic fox] \
Script: [CanVulpAncestralFasta.sh]()

### STEP1: Find a ‘global estimate’ of the SFS
#### STEP1a: Estimate the site allele frequency likelihood
*Note: If you don't have the ancestral state, you can instead estimate the folded SFS. This is done by supplying the -anc with the reference genome and applying -fold 1 to realSFS. (aka do not supply a -ref)*

[/home/sophiepq/GrayFox/angsd/SFS/EastFullK3HC_doSaf.sh](https://github.com/squisquater/Gray-Fox-Hybrid-Zone-Dynamics/blob/main/Demographic-History/SFS/EastFullK3HC_doSaf.sh) \
[/home/sophiepq/GrayFox/angsd/SFS/TexasFullK3HC_doSaf.sh](https://github.com/squisquater/Gray-Fox-Hybrid-Zone-Dynamics/blob/main/Demographic-History/SFS/TexasFullK3HC_doSaf.sh) \
[/home/sophiepq/GrayFox/angsd/SFS/EastPureK3HC_doSaf.sh](https://github.com/squisquater/Gray-Fox-Hybrid-Zone-Dynamics/blob/main/Demographic-History/SFS/EastPureK3HC_doSaf.sh) \
[/home/sophiepq/GrayFox/angsd/SFS/TexasPureK3HC_doSaf.sh](https://github.com/squisquater/Gray-Fox-Hybrid-Zone-Dynamics/blob/main/Demographic-History/SFS/TexasPureK3HC_doSaf.sh)

#### STEP1b: Generate a 1-D site frequency spectrum --> [1D-SFS.sh](https://github.com/squisquater/Gray-Fox-Hybrid-Zone-Dynamics/blob/main/Demographic-History/SFS/1D-SFS.sh)

## Generate 2-D Site FrequencySpectrum 
* [2D-SFS-ET-FullK3HC.sh](https://github.com/squisquater/Gray-Fox-Hybrid-Zone-Dynamics/blob/main/Demographic-History/SFS/2D-SFS-ET-FullK3HC.sh)
* [2D-SFS-ET-PureK3HC.sh](https://github.com/squisquater/Gray-Fox-Hybrid-Zone-Dynamics/blob/main/Demographic-History/SFS/2D-SFS-ET-PureK3HC.sh)

Need to run some brief R code to convert these into 2D matrices.

```
sfsETpure<-scan("2D-SFS-ET-PureK3HC.sfs")
sfsETfull<-scan("2D-SFS-ET-FullK3HC.sfs")

nEfull=69
nTfull=87
nEpure=61 
nTpure=27

SFSmatrixETfull.byrow <- matrix(sfsETfull,nrow=((2*nEfull)+1),ncol=((2*nTfull)+1), byrow=T)
SFSmatrixETpure.byrow <- matrix(sfsETpure,nrow=((2*nEpure)+1),ncol=((2*nTpure)+1), byrow=T)

SFSmatrixETfull.byrow <- round(SFSmatrixETfull.byrow,0)
SFSmatrixETpure.byrow <- round(SFSmatrixETpure.byrow,0)

rownames(SFSmatrixETfull.byrow) <- rownames(SFSmatrixETfull.byrow, do.NULL = FALSE, prefix = "d_")
colnames(SFSmatrixETfull.byrow) <- colnames(SFSmatrixETfull.byrow, do.NULL = FALSE, prefix = "d_")
rownames(SFSmatrixETpure.byrow) <- rownames(SFSmatrixETpure.byrow, do.NULL = FALSE, prefix = "d_")
colnames(SFSmatrixETpure.byrow) <- colnames(SFSmatrixETpure.byrow, do.NULL = FALSE, prefix = "d_")

write.table(SFSmatrixETfull.byrow, file = "sfsETfull_jointDAFpop1_0.obs", sep="\t", row.names=T, col.names=T, quote=F)
write.table(SFSmatrixETpure.byrow, file = "sfsETpure_jointDAFpop1_0.obs", sep="\t", row.names=T, col.names=T, quote=F)

#Note I still needed to make a couple manual tweaks to this input file. I added 
#"1 observation" to the first line
#and added a tab indentation to the column header line (line2)
```
## Run Fastsimcoal2 Models


### ET Pure - no migration
sfsETpure-nomig_jointDAFpop1_0.obs
sfsETpure-nomig.tpl
```
//Parameters for the coalescence simulation program : simcoal.exe
2 samples to simulate :
//Population effective sizes (number of genes)
NPOP1
NPOP2
//Samples sizes and samples age
122
54
//Growth rates: negative growth implies population expansion
0
0
//Number of migration matrices : 0 implies no migration between demes
0
//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index
1 historical event
TDIV 0 1 1 RESIZE 0 1
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per gen recomb and mut rates
FREQ 1 0 4.5e-9 OUTEXP
```

sfsETpure-nomig.est
```
// Priors and rules file
// *********************

[PARAMETERS]
//#isInt? #name   #dist.#min  #max
//all N are in number of haploid individuals
1  ANCSIZE     unif     100  5e5   output
1  NPOP1       unif     100  5e5   output
1  NPOP2       unif     100  5e5   output
0  N1M21       logunif  1e-2 20       hide
0  N2M12       logunif  1e-2 20       hide
1  TDIV        unif     1000   1e6   output

[RULES]

[COMPLEX PARAMETERS]
0  RESIZE = ANCSIZE/NPOP2     hide
```
### ET Pure - constant migration
sfsETpure-constantmig_jointDAFpop1_0.obs
sfsETpure-constantmig.tpl
```
//Parameters for the coalescence simulation program : simcoal.exe
2 samples to simulate :
//Population effective sizes (number of genes)
NPOP1
NPOP2
//Samples sizes and samples age
122
54
//Growth rates: negative growth implies population expansion
0
0
//Number of migration matrices : 0 implies no migration between demes
2
//Migration matrix 0
0 MIG21
MIG12 0
//Migration matrix 1
0 0
0 0
//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index
1 historical event
TDIV 0 1 1 RESIZE 0 1
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per gen recomb and mut rates
FREQ 1 0 4.5e-9 OUTEXP
```

sfsETpure-constantmig.est
```
// Priors and rules file
// *********************

[PARAMETERS]
//#isInt? #name   #dist.#min  #max
//all N are in number of haploid individuals
1  ANCSIZE     unif     100  5e5   output
1  NPOP1       unif     100  5e5   output
1  NPOP2       unif     100  5e5   output
0  N1M21       logunif  1e-2 20       hide
0  N2M12       logunif  1e-2 20       hide
1  TDIV        unif     1000   1e6   output

[RULES]

[COMPLEX PARAMETERS]
0  RESIZE = ANCSIZE/NPOP2     hide
0  MIG12  = N1M21/NPOP1       output
0  MIG21  = N2M12/NPOP2       output
```
### ET Pure - recent gene flow
sfsETpure-migstop_jointDAFpop1_0.obs
sfsETpure-migstop.tpl
```
//Parameters for the coalescence simulation program : simcoal.exe
2 samples to simulate :
//Population effective sizes (number of genes)
NPOP1
NPOP2
//Samples sizes and samples age
122
54
//Growth rates: negative growth implies population expansion
0
0
//Number of migration matrices : 0 implies no migration between demes
2
//Migration matrix 0
0 MIG21
MIG12 0
//Migration matrix 1
0 0
0 0
//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index
2 historical event
TMIG 0 0 0 1 0 1
TDIV 0 1 1 RESIZE 0 1
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per gen recomb and mut rates
FREQ 1 0 4.5e-9 OUTEXP
```

sfsETpure-migstop.est
```
// Priors and rules file
// *********************

[PARAMETERS]
//#isInt? #name   #dist.#min  #max
//all N are in number of haploid individuals
1  ANCSIZE     unif     100  5e5   output
1  NPOP1       unif     100  5e5   output
1  NPOP2       unif     100  5e5   output
0  N1M21       logunif  1e-2 20       hide
0  N2M12       logunif  1e-2 20       hide
1  TDIV        unif     1000   1e6   output
1  TMIG        unif     10   1e6   output       

[RULES]
TDIV>TMIG

[COMPLEX PARAMETERS]
0  RESIZE = ANCSIZE/NPOP2     hide
0  MIG12  = N1M21/NPOP1       output
0  MIG21  = N2M12/NPOP2       output
```


## Run fastsimcoal for all analyses
*sfsETpure-nomig.sh*
```
#!/usr/bin/env bash
#SBATCH --job-name=sfsETpure-nomig
#SBATCH --array=1-100%33
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --time 10:00:00
#SBATCH --mem=2GB
#SBATCH -o /home/sophiepq/GrayFox/angsd/SFS/fastsimcoal/sfsETpure-nomig/sfsETpure-nomig_%A_%a.out
#SBATCH -e /home/sophiepq/GrayFox/angsd/SFS/fastsimcoal/sfsETpure-nomig/sfsETpure-nomig_%A_%a.out
#SBATCH -p high

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
STARTTIME=$(date +"%s")

MODEL=sfsETpure-nomig
RUN=$(sed "${SLURM_ARRAY_TASK_ID}q;d" /home/sophiepq/GrayFox/angsd/SFS/fastsimcoal/sfsETpure-nomig/100runs.txt | cut -f1)

echo "running run $RUN of the fastsimcoal analysis on $MODEL model"

cd /home/sophiepq/GrayFox/angsd/SFS/fastsimcoal/$MODEL
mkdir Run$RUN
cp ./$MODEL.tpl ./Run$RUN
cp ./$MODEL.est ./Run$RUN
cp ./${MODEL}_jointDAFpop1_0.obs ./Run$RUN
cd ./Run$RUN

~/bin/fsc26_linux64/fsc26 -t $MODEL.tpl -n200000 -d -e $MODEL.est -M -L 75 -c6 -q

ENDTIME=$(date +%s)
TIMESPEND=$(($ENDTIME - $STARTTIME))
((sec=TIMESPEND%60,TIMESPEND/=60, min=TIMESPEND%60, hrs=TIMESPEND/60))
timestamp=$(printf "%d:%02d:%02d" $hrs $min $sec)
echo "Took $timestamp hours:minutes:seconds to complete..."
```


*sfsETpure-constantmig.sh*
```
#!/usr/bin/env bash
#SBATCH --job-name=sfsETpure-constantmig
#SBATCH --array=1-100%33
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --time 20:00:00
#SBATCH --mem=2GB
#SBATCH -o /home/sophiepq/GrayFox/angsd/SFS/fastsimcoal/sfsETpure-constantmig/sfsETpure-constantmig_%A_%a.out
#SBATCH -e /home/sophiepq/GrayFox/angsd/SFS/fastsimcoal/sfsETpure-constantmig/sfsETpure-constantmig_%A_%a.out
#SBATCH -p high

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
STARTTIME=$(date +"%s")

MODEL=sfsETpure-constantmig
RUN=$(sed "${SLURM_ARRAY_TASK_ID}q;d" /home/sophiepq/GrayFox/angsd/SFS/fastsimcoal/sfsETpure-constantmig/100runs.txt | cut -f1)

echo "running run $RUN of the fastsimcoal analysis on $MODEL model"

cd /home/sophiepq/GrayFox/angsd/SFS/fastsimcoal/$MODEL
mkdir Run$RUN
cp ./$MODEL.tpl ./Run$RUN
cp ./$MODEL.est ./Run$RUN
cp ./${MODEL}_jointDAFpop1_0.obs ./Run$RUN
cd ./Run$RUN

~/bin/fsc26_linux64/fsc26 -t $MODEL.tpl -n200000 -d -e $MODEL.est -M -L 75 -c6 -q

ENDTIME=$(date +%s)
TIMESPEND=$(($ENDTIME - $STARTTIME))
((sec=TIMESPEND%60,TIMESPEND/=60, min=TIMESPEND%60, hrs=TIMESPEND/60))
timestamp=$(printf "%d:%02d:%02d" $hrs $min $sec)
echo "Took $timestamp hours:minutes:seconds to complete..."
```

*sfsETpure-migstop.redo.sh*
```
#!/usr/bin/env bash
#SBATCH --job-name=sfsETpure-migstop.redo
#SBATCH --array=1-100%33
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --time 20:00:00
#SBATCH --mem=2GB
#SBATCH -o /home/sophiepq/GrayFox/angsd/SFS/fastsimcoal/sfsETpure-migstop.redo/sfsETpure-migstop.redo_%A_%a.out
#SBATCH -e /home/sophiepq/GrayFox/angsd/SFS/fastsimcoal/sfsETpure-migstop.redo/sfsETpure-migstop.redo_%A_%a.out
#SBATCH -p high

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
STARTTIME=$(date +"%s")

MODEL=sfsETpure-migstop.redo
RUN=$(sed "${SLURM_ARRAY_TASK_ID}q;d" /home/sophiepq/GrayFox/angsd/SFS/fastsimcoal/sfsETpure-migstop.redo/100runs.txt | cut -f1)

echo "running run $RUN of the fastsimcoal analysis on $MODEL model"

cd /home/sophiepq/GrayFox/angsd/SFS/fastsimcoal/$MODEL
mkdir Run$RUN
cp ./$MODEL.tpl ./Run$RUN
cp ./$MODEL.est ./Run$RUN
cp ./${MODEL}_jointDAFpop1_0.obs ./Run$RUN
cd ./Run$RUN

~/bin/fsc26_linux64/fsc26 -t $MODEL.tpl -n200000 -d -e $MODEL.est -M -L 75 -c6 -q

ENDTIME=$(date +%s)
TIMESPEND=$(($ENDTIME - $STARTTIME))
((sec=TIMESPEND%60,TIMESPEND/=60, min=TIMESPEND%60, hrs=TIMESPEND/60))
timestamp=$(printf "%d:%02d:%02d" $hrs $min $sec)
echo "Took $timestamp hours:minutes:seconds to complete..."
```

## Extract the highest likelihood model from each of the 100 replicate runs

To find the best run, i.e. the run with the highest likelihood, or better the smallest difference between the maximum possible likelihood (MaxObsLhood) and the obtained likelihood (MaxEstLhood), we can check the .bestlhoods files.

```
cat Run{1..91}/sfsETpure-nomig/sfsETpure-nomig.bestlhoods | grep -v MaxObsLhood | awk '{print NR,$8}' | sort -k 2

cat Run{1..100}/sfsETAll-nomig/sfsETAll-nomig.bestlhoods | grep -v MaxObsLhood

cat Run{1..40}/sfsETpure-nomig/sfsETpure-nomig.bestlhoods | grep -v MaxObsLhood

cat Run{1..91}/${PREFIX}/${PREFIX}.bestlhoods | grep -v MaxObsLhood | awk '{print NR,$8}' | sort -k 2


```

```
#!/bin/bash

# Bash script by Joana Meier to select the best fastsimcoal run of multiple runs under the same demographic model
# The script expects output files of different runs to be found in folders starting with run 

m=-1000000000000000
p=$1
c=0
best="xxxx"
diff=0
diffBest=0

for i in run*;
do
 a=$(ls | grep ${p}".tpl" | sed s'/.tpl//')

 # if the file is in the run directory
 if [ -e $i/$a.bestlhoods ]
 then
   l=$(cat $i/$a.bestlhoods | awk '{print $(NF-1)}' | tail -1 | cut -f 1 -d ".")
   diff=$(cat $i/$a.bestlhoods | awk '{print $NF-$(NF-1)}' | tail -1 | cut -f 1 -d ".")
   ((c++))
 else
   # if the file is in a subdirectory
   if [ -e $i/$a/$a.bestlhoods ]
   then
     l=$(cat $i/$a/$a.bestlhoods | awk '{print $(NF-1)}' | tail -1 | cut -f 1 -d ".")
     diff=$(cat $i/$a/$a.bestlhoods | awk '{print $NF-$(NF-1)}' | tail -1 | cut -f 1 -d ".")
     ((c++))
   else
     echo "no .bestlhoods file found in "$i
   fi

 fi


 if [ $l -gt $m ]
 then
   m=$l; x=$i; best=$i;
   diffBest=$diff
 fi
done

if [ -z ${x+1} ]
then echo "Error: No run with lik>-1000000000000000"
else mkdir bestrun
cp $x/* ./bestrun/
cp $x/$a/* ./bestrun/
fi

echo -e "\n"$c" bestlhoods files found, "$best" with "$diffBest" Lhood diff fits best."
```



##### IGNORE BELOW FOR THE MOMENT 
#### Several Other Datasets That Could be Used

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

[/home/sophiepq/GrayFox/angsd/SFS/EastFullK3HC_doSaf.sh]() **running** \
[/home/sophiepq/GrayFox/angsd/SFS/WestFullK3HC_doSaf.sh]() **running** \
[/home/sophiepq/GrayFox/angsd/SFS/TexasFullK3HC_doSaf.sh]() **running** \
[/home/sophiepq/GrayFox/angsd/SFS/CaliforniaFullK3HC_doSaf.sh]() **running** 

[/home/sophiepq/GrayFox/angsd/SFS/EastPureAll_doSaf.sh]() **to run** \
[/home/sophiepq/GrayFox/angsd/SFS/WestPureAll_doSaf.sh]() **to run** \
[/home/sophiepq/GrayFox/angsd/SFS/TexasK3PureAll_doSaf.sh]() **to run** \
[/home/sophiepq/GrayFox/angsd/SFS/AllCaliforniaPureAll_doSaf.sh]() **to run** \
[/home/sophiepq/GrayFox/angsd/SFS/N.CAPureAll_doSaf.sh]() **to run** 


#### STEP1b: Generate a 1-D site frequency spectrum --> [1D-SFS.sh](https://github.com/squisquater/Gray-Fox-Hybrid-Zone-Dynamics/blob/main/Demographic-History/SFS/1D-SFS.sh)

## Generate 2-D Site FrequencySpectrum 
* [2D-SFS-EW.sh]()
* [2D-SFS-ET.sh]()
* [2D-SFS-EC.sh]()
* [2D-SFS-CT.sh]()
* [2D-SFS-ET-PureK3HC.sh](https://github.com/squisquater/Gray-Fox-Hybrid-Zone-Dynamics/blob/main/Demographic-History/SFS/2D-SFS-ET-PureK3HC.sh)

## Generate 3-D Site FrequencySpectrum 
* [3D-SFS-ETC.sh]()



