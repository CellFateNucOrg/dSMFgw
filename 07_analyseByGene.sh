#! /bin/bash
#SBATCH --mail-user=jennifer.semple@izb.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="EMs_euc"
#SBATCH --partition=epyc2
#SBATCH --mem-per-cpu=4G
#SBATCH --time=0-02:00:00
#SBATCH --output=slurm-%x-%A_%a.out
#SBATCH --error=slurm-%x-%A_%a.out 

module add vital-it
module load R/3.6.1


# Call methylation and do plots
source ./varSettings.sh

Rscript --no-save --no-restore --verbose plotByGene.R 2>&1

