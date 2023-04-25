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
For the PSMC plotting code, it is important to include the flag -R, which will give you the data (population size and years) to plot in ggplot2. Without this flag, PSMC will just output a figure that it will make, which is hard to adjust and plot multiple individual's demographic history trajectories on.

Generate a conda environment which has necessary packages to run the plot.py script.
```
conda create --name PSMC
conda install -c bioconda gnuplot
conda activate PSMC
```
