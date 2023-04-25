### PSMC Analysis
**STEP0:** Subset autosomes [[subset_autosomes.sh]](https://github.com/squisquater/Gray-Fox-Hybrid-Zone-Dynamics/blob/main/Demographic-History/PSMC/subset_autosomes.sh) \
**STEP1:** Convert bam to fastq file [[bam2fastq.sh]](https://github.com/squisquater/Gray-Fox-Hybrid-Zone-Dynamics/blob/main/Demographic-History/PSMC/bam2fastq.sh) \
**STEP2:** Convert fastq to psmc fasta file [[fastq2psmcfa.sh]](https://github.com/squisquater/Gray-Fox-Hybrid-Zone-Dynamics/blob/main/Demographic-History/PSMC/fastq2psmcfa.sh) \
**STEP3:** Run PSMC for all files using Wang et al (2020) parameters [[runPSMC.sh]]() \
**STEP4:** Bootstrapping 
> Split fasta files for bootstrapping [[bootstrap_split.sh]](https://github.com/squisquater/Gray-Fox-Hybrid-Zone-Dynamics/blob/main/Demographic-History/PSMC/bootstrap_split.sh) \
> Run bootstrapping [[bootstrap.sh]](https://github.com/squisquater/Gray-Fox-Hybrid-Zone-Dynamics/blob/main/Demographic-History/PSMC/bootstrap.sh) 

**STEP5:** Merge all .psmc files together [[mergePSMC.sh]](https://github.com/squisquater/Gray-Fox-Hybrid-Zone-Dynamics/blob/main/Demographic-History/PSMC/mergePSMC.sh) \
**STEP6:** Plot PSMC trajectories [[plotPSMC.sh]](https://github.com/squisquater/Gray-Fox-Hybrid-Zone-Dynamics/blob/main/Demographic-History/PSMC/plotPSMC.sh)
> For the PSMC plotting code, it is important to include the flag -R, which will give you the data (population size and years) to plot in ggplot2. 
**STEP7:** Combine all the bootstrapped replicates so I can plot them together [[mergePSMCplots.sh]]()
