#! /bin/bash
##SBATCH --mail-user=bolaji.isiaka@izb.unibe.ch
##SBATCH --mail-type=end,fail
#SBATCH --job-name="dSMF_plots"
#SBATCH --time=2-00:00:00
#SBATCH --partition=all
#SBATCH --mem-per-cpu=12G

##SBATCH --tmp=64G

#echo $TMPDIR "is temp dir"

module add vital-it
module load R/3.6.1
module add UHTS/Analysis/MultiQC/1.7;
module add UHTS/Analysis/samtools/1.8;

# Collect various QC data produced by previous script together
#multiqc -f ./qc

# Call methylation and do plots
source ./varSettings.sh


#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#Rscript dSMFseqAnalysis_singleMolecule.R $SLURM_ARRAY_TASK_ID $SLURM_CPUS_PER_TASK
Rscript dSMFseqAnalysis_metagene.R

