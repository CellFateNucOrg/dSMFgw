#! /usr/bin/bash

## Allocate resources
#SBATCH --time=0-24:00:00
#SBATCH --array=1-2

#SBATCH --mail-user=jennifer.semple@izb.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="bwaMeth_dSMF"
#SBATCH --cpus-per-task=4
#SBATCH --partition=all
#SBATCH --mem-per-cpu=8G

module add vital-it
module load R/3.5.1
module add UHTS/Quality_control/fastqc/0.11.5      #fastqc
module add UHTS/Quality_control/cutadapt/1.13     #cutadapt
module add UHTS/Analysis/trimmomatic/0.36;
module add UHTS/Analysis/samtools/1.8;
module add UHTS/Analysis/picard-tools/2.18.11;
module add UHTS/Quality_control/qualimap/2.2.1;
module add UHTS/Aligner/bwa/0.7.17;
#source activate bwaMeth

# read in the run specific settings
source ./varSettings.sh

# get index for list of samples
let i=$SLURM_ARRAY_TASK_ID-1

# do QC and map 
./alnBwaMeth.sh ${sampleNames[$i]} $SLURM_CPUS_PER_TASK

