#! /usr/bin/bash

## Allocate resources
#SBATCH --time=3-00:00:00
#SBATCH --array=1-2
#SBATCH --mail-user=jennifer.semple@izb.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="bismark_dSMF"
#SBATCH --cpus-per-task=4
#SBATCH --partition=all
#SBATCH --mem-per-cpu=8G
#SBATCH --tmp=32G

module add vital-it;
module load R/3.5.1;
module add UHTS/Quality_control/fastqc/0.11.5;      #fastqc
module add UHTS/Quality_control/cutadapt/1.13;     #cutadapt
module add UHTS/Analysis/trimmomatic/0.36;
module add UHTS/Analysis/samtools/1.8;
module add UHTS/Analysis/picard-tools/2.18.11;
module add UHTS/Quality_control/qualimap/2.2.1;
module add UHTS/Aligner/bowtie2/2.3.4.1;
#module add UHTS/Aligner/bwa/0.7.17;
#module add UHTS/Analysis/BBMap/37.82;
module add UHTS/Analysis/bamtools/2.4.1;
#module add UHTS/Analysis/MultiQC/1.7;
#source activate bwaMeth

# read in the run specific settings
source ./varSettings.sh

# get index for list of samples
let i=$SLURM_ARRAY_TASK_ID-1

# do QC and map 
if [[ "$dataType" == gw ]]
then
	./alnBwaMeth_gw.sh ${sampleNames[$i]} ${testGroups[$i]} ${SLURM_CPUS_PER_TASK}
        echo "processing genome wide library"
elif [[ "$dataType" == amp ]]
then
	./alnBwaMeth_amp.sh ${sampleNames[$i]} ${testGroups[$i]} ${SLURM_CPUS_PER_TASK}
        echo "processing amplicon library"
else
        echo "ERROR: dataType in varSettings.sh must be either "amp" (amplicon library) or "gw" (genome-wide library)"
fi
