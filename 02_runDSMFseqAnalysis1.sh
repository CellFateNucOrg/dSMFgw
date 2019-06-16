#! /bin/bash
#SBATCH --mail-user=jennifer.semple@izb.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="feb18"
#SBATCH --time=0-00:10:00
#SBATCH --partition=debug
#SBATCH --mem-per-cpu=64G

module add vital-it
module load R/3.5.1
module add UHTS/Analysis/MultiQC/1.7;
module add UHTS/Analysis/samtools/1.8;

# Collect various QC data produced by previous script together
#multiqc -f ./qc

# Call methylation and do plots
source ./varSettings.sh

Rscript dSMFseqAnalysis.R 

