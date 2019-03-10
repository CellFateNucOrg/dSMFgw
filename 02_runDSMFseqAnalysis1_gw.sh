#! /bin/bash
#SBATCH --mail-user=jennifer.semple@izb.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="gw_dSMF"
#SBATCH --time=1-12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --partition=all
#SBATCH --mem-per-cpu=8G

module add vital-it
module load R/3.5.1

source ./varSettings.sh

Rscript dSMFseqAnalysis1_gw_aln.R ${genomefile}

