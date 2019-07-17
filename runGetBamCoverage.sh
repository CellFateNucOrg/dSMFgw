#! /usr/bin/bash

## Allocate resources
#SBATCH --time=0-02:00:00
#SBATCH --mail-user=jennifer.semple@izb.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="getCov"
#SBATCH --cpus-per-task=4
#SBATCH --partition=all
#SBATCH --mem-per-cpu=8G


module load vital-it;
module add UHTS/Analysis/deepTools/2.5.4;

./getBamCoverage.sh aln/*.noOL.bam
