#! /usr/bin/bash

## Allocate resources
#SBATCH --time=0-00:10:00
#SBATCH --mail-user=jennifer.semple@izb.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="formatUCSC"
#SBATCH --cpus-per-task=4
#SBATCH --partition=all
#SBATCH --mem-per-cpu=8G



./formatForUCSC.sh gwDSMF16-20_Feb "Genome wide dSMF dS016 & dS020, Feb-2019" ./bigwig/*.bw

