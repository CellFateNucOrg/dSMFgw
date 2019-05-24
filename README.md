# dSMFgw

Pipeline for analysing genome-wide dSMF data. This pipeline aligns dSMF reads with bwa-meth rather than QuasR, as the quality and yield of the alignments are much better. Only the bam files are imported into QuasR for methylation calling.

## Installing bwa-meth

```
# installation instructions for bwa-meth. based on https://github.com/brentp/bwa-meth
conda create --name bwameth python=3.7

source activate bwameth

pip install toolshed

wget https://github.com/brentp/bwa-meth/archive/master.zip
unzip master.zip
cd bwa-meth-master/

# local install without sudo privilages
python setup.py install --user

# bwa-meth will be at: ~/.local/bin/bwameth.py
# can run using

conda activate bwameth
python bwameth.py


# Note: you must have samtools and bwa-mem modules activated
module add vital-it
module add UHTS/Analysis/samtools/1.8
source activate bwameth

#Before using it you must index the genomecon


```

## Installing MethylDackel
```

conda activate bwameth
conda install -c bioconda methyldackel
conda deactivatemeth
````

Note, this is only version 3.0.
Compiling the dev version didn't seem to work but didn't give an error either:
```
module load /software/UHTS/Analysis/HTSlib/1.9;

make install CFLAGS="-O3 -Wall -I/software/UHTS/Analysis/HTSlib/1.9/include" LIBS="-L/software/UHTS/Analysis/HTSlib/1.9/lib" prefix=/home/ubelix/izb/semple/mySoftware/methyldackel

