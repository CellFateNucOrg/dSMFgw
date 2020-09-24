#! /bin/bash
#SBATCH --mail-user=jennifer.semple@izb.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="EMmulti_euclid"
#SBATCH --array=2-3
#SBATCH --cpus-per-task=4
#SBATCH --partition=all
#SBATCH --mem-per-cpu=8G
#SBATCH --time=3-00:00:00
#SBATCH --output=slurm-%A_%a-%x.out
#SBATCH --error=slurm-%A_%a-%x.out 

module add vital-it
module load R/3.6.1


# Call methylation and do plots
source ./varSettings.sh

echo "task id is: " $SLURM_ARRAY_TASK_ID
Rscript --no-save --no-restore --verbose clusterMultipleGenes.R $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_MAX $SLURM_CPUS_PER_TASK 2>&1

