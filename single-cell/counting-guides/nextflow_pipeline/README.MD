
# Runing our starSolo nextflow module to count guides
This workflow will receive the whitelist, fastq and metadata and generate as output an AnnData (cells x guides) file using starSolo.

```bash
#iNSTALL CONDA 
conda update -n base conda -y 
conda install -n base conda-libmamba-solver -y
conda config --set solver libmamba #this make conda very fast
#install nextflow
conda install bioconda::nextflow

```

```bash

#unzip the  data sample
tar -xvzf example_data/data.tar.gz
#This data is very small should be doable running in your desktop or small server
```

```bash

#Edit the input.config config file to match the directory that you extracted your data.
#The paths should correspond to your computer/server paths

```
```bash
#Execute the nextflow from the project directory

nextflow main.nf -config input.config -resume


```

## Directory structure:

### bin 

 all custom script should be distributed inside this directory. Make sure allow  for execution of your  after a script creation ex: __chmod +x script.sh__

### conda_envs 

 The conda file with your enviroment installation commands image. If you need a docker environment create a docker_envs directory

### processes
Directory with the nextflow processes.

### example_data
A directory with your data. The input.config file should point to the pathway of these files

### input.config
 Document with the parameters being passed as input for your nextflow module

### main.nf
Your module workflow






