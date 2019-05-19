#! /bin/bash
#SBATCH --mail-user=jennifer.semple@izb.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="dSMF_Ccall"
#SBATCH --time=2-00:00:00
##SBATCH --cpus-per-task=4
#SBATCH --partition=all
#SBATCH --mem-per-cpu=64G
##SBATCH --array=1-2

module add vital-it
module load R/3.5.1
module add UHTS/Analysis/MultiQC/1.7;
module add UHTS/Analysis/samtools/1.8;

# Collect various QC data produced by previous script together
#multiqc ./fastQC

# Call methylation and do plots
source ./varSettings.sh

Rscript dSMFseqAnalysis.R 

