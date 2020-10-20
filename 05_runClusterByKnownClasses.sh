#! /bin/bash
#SBATCH --mail-user=jennifer.semple@izb.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="EMm_euc"
#SBATCH --array=1-3
#SBATCH --cpus-per-task=4
#SBATCH --partition=all
#SBATCH --mem-per-cpu=12G
#SBATCH --time=0-01:00:00
#SBATCH --output=slurm-%x-%A_%a.out
#SBATCH --error=slurm-%x-%A_%a.out 

module add vital-it
module load R/3.6.1


# Call methylation and do plots
source ./varSettings.sh

echo "task id is: " $SLURM_ARRAY_TASK_ID
Rscript --no-save --no-restore --verbose clusterByKnownClasses_multigene.R $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_MAX $SLURM_CPUS_PER_TASK 2>&1

