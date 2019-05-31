# dSMFgw

Pipeline for analysing genome-wide dSMF data. This pipeline aligns dSMF reads with bwa-meth rather than QuasR, as the quality and yield of the alignments are much better. Only the bam files are imported into QuasR for methylation calling.

## Installing bwa-meth

```
# installation instructions for bwa-meth. based on https://github.com/brentp/bwa-meth
conda create --name bwaMeth python=3.7

source activate bwaMeth

pip install toolshed

wget https://github.com/brentp/bwa-meth/archive/master.zip
unzip master.zip
cd bwa-meth-master/

# local install without sudo privilages
python setup.py install --user

# bwa-meth will be at: ~/.local/bin/bwameth.py
# add a line to .bashrc to point to it:
cp ~/.bashrc ~/.bashrc_backup
echo "export BWAMETH=${HOME}/.local/bin/bwameth.py" >> ~/.bashrc

# to be able to activate environments from inside a slurm script you need to add
# the path to the activate script to .bashrc:
export CONDA_ACTIVATE=/home/ubelix/izb/semple/anaconda3/bin/activate
# then in the script use
source $CONDA_ACTIVATE
conda activate bwameth

# Note: you must have samtools and bwa-mem modules activated

#Before using it you must index the genome


#leave environment
conda deactivate
```

## Installing MethylDackel
```
#conda create --name methyldackel python=3.7
conda activate bwameth
conda install -c bioconda methyldackel
conda deactivate
````

Note, this is only version 3.0.
Compiling the dev version didn't seem to work but didn't give an error either:
```
module load /software/UHTS/Analysis/HTSlib/1.9;

make install CFLAGS="-O3 -Wall -I/software/UHTS/Analysis/HTSlib/1.9/include" LIBS="-L/software/UHTS/Analysis/HTSlib/1.9/lib" prefix=/home/ubelix/izb/semple/mySoftware/methyldackel

