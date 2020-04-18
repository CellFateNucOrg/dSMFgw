#! /bin/bash
#SBATCH --mail-user=jennifer.semple@izb.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="dSMF_EMclust"
#SBATCH --array=1-8
#SBATCH --cpus-per-task=4
#SBATCH --partition=all
#SBATCH --mem-per-cpu=4G

module add vital-it
module load R/3.6.1


# Call methylation and do plots
source ./varSettings.sh

Rscript --no-save --no-restore --verbose clusterSingleGenes.R $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_MAX $SLURM_CPUS_PER_TASK 2>&1

