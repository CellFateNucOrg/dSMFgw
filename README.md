# dSMFgw

Pipeline for analysing genome-wide dSMF data. This pipeline aligns dSMF reads with bwa-meth rather than QuasR, as the quality and yield of the alignments are much better. Only the bam files are imported into QuasR for methylation calling.

## Installing bwa-meth

```
# installation instructions for bwa-meth. based on https://github.com/brentp/bwa-meth
conda create --name bwaMeth python=2.7

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

# Note: you must have samtools and bwa-mem modules activated

#Before using it you must index the genome


#leave environment
source deactivate bwameth
```

## Installing MethylDackel
conda create --name methyldackel python=3.7
conda activate methyldackel
conda install -c bioconda methyldackel
conda deactivate


