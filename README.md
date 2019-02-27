# bwaMeth

Aligning dSMF reads with bwa-meth rather than QuasR.

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

