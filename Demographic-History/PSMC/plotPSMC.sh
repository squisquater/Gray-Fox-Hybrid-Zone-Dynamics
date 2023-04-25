!/bin/bash -l
#SBATCH --job-name=PSMC_S19-3133_plot
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --time 10:00
#SBATCH --mem=100M
#SBATCH -p bmm
#SBATCH -A ctbrowngrp
#SBATCH -o PSMC_S19-3133_plot.out
#SBATCH -e PSMC_S19-3133_plot.err


/group/ctbrowngrp2/hennelly/hennelly/bin/psmc/utils/psmc_plot.pl -R -u 4.5e-09 -g 2 S19_3133 S19_3133_combined.psmc


## Note this will just generate output .txt files for each replicate so you can plot them in ggplot. This will not generate plots for you as is done in some PSMC tutorials
