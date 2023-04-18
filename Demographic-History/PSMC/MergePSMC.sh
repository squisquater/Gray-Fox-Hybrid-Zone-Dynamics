#This is super quick so I just ran this interactively.
#You're essentially combining the initial (unsplit) PSMC analysis with all the bootstrapped replicates.

srun --partition=high --time=2:00:00 --ntasks=1 --cpus-per-task=4 --mem=1G --pty /bin/bash -l

cd /group/ctbrowngrp2/sophiepq/GrayFoxWGS/GrayFox1/angsd/PSMC/PSMC_files_new/bootstrap/PSMC_S19_6655
cat S19_6655_GFaddon_autosomes.psmc run-*.psmc > S19_6655_combined.psmc

cd /group/ctbrowngrp2/sophiepq/GrayFoxWGS/GrayFox1/angsd/PSMC/PSMC_files_new/bootstrap/PSMC_S19_3133
cat S19_3133_GFaddon_autosomes.psmc run-*.psmc > S19_3133_combined.psmc
