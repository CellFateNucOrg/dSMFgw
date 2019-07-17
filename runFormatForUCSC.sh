#! /usr/bin/bash

## Allocate resources
#SBATCH --time=0-00:10:00
##SBATCH --mail-user=jennifer.semple@izb.unibe.ch
##SBATCH --mail-type=end,fail
#SBATCH --job-name="formatUCSC"
#SBATCH --cpus-per-task=1
#SBATCH --partition=all
#SBATCH --mem-per-cpu=2G



./formatForUCSC.sh gwDSMFv16v20_Comb "Genome wide dSMF dS016 & dS020, Feb & May 2019 combined" ./bigwig/*.bw

