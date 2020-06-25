# dSMFgw

Pipeline for analysing genome-wide dSMF data. This pipeline aligns dSMF reads with bwa-meth rather than QuasR, as the quality and yield of the alignments are much better. Only the bam files are imported into QuasR for methylation calling.

## Installing bwa-meth

```
# installation instructions for bwa-meth. based on https://github.com/brentp/bwa-meth
conda create --name bwaMeth python=3.7

conda activate bwaMeth

pip install toolshed
conda install -c bioconda samtools

wget https://github.com/brentp/bwa-meth/archive/master.zip
unzip master.zip
cd bwa-meth-master/

# local install without sudo privilages
python setup.py install --user

# bwa-meth will be at: ~/.local/bin/bwameth.py
# add a line to .bashrc to point to it:
cp ~/.bashrc ~/.bashrc_backup
echo "export BWAMETH=${HOME}/.local/bin/bwameth.py" >> ~/.bashrc

# Note: you must have bwa-mem module activated or install it separately

#leave environment
conda deactivate
```

## Installing MethylDackel
```
conda activate bwaMeth
conda install -c bioconda methyldackel
````

## Install bamUtil
See info at: https://genome.sph.umich.edu/wiki/BamUtil#Building
https://github.com/statgen/bamUtil
```
wget https://github.com/statgen/bamUtil/archive/master.tar.gz
```
Unzip it with tar -xzvf, then install in your local software directory (e.g.${HOME}/mySoftware/):
```
make cloneLib
make
make install INSTALLDIR=${HOME}/mySoftware/bamUtil
```
You will need to add the path of the executable to the varSettings.sh file

## Install R libraries
Check the code at the beginning of 02_ script to see which libraries to install.
You need to install the methMatrix library from github as follows:
```
devtools::install_github("jsemple19/methMatrix")


## varSettings.sh
most of what you have to change goes into the varSettings.sh file. An example file comes with the repository, copy it and change the name:
```
cp varSettings_example.sh varSettings.sh
```
Then use a text editor to change verSettings.sh

The only other changes you will have to make is to the number of array jobs in 01_ script: The number of jobs should be the same as the number of libraries you sent to sequence.

# Running the scripts
run the scripts using that sbatch wrapper scripts 01_, 02a_, 02b_ and 02c_ in the order indicated by the number
