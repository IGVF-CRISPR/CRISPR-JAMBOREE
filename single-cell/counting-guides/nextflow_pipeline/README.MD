
# Runing  starSolo nextflow module to count guides
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
cd CRISPR-JAMBOREE/single-cell/counting-guides/nextflow_pipeline/example_data  
tar -xvzf -k data.tar.gz
#This data is very small should be doable running in your desktop or small server
```

```bash

#make the bin files executable
chmod +x bin/*    
```


```bash

#Edit the input.config config file to match the directory that you extracted your data.
#The paths should correspond to your computer/server paths

```
```bash
#Execute the nextflow from the project directory

nextflow main.nf -config input.config -resume -with-conda


```
## Inputs
- GUIDES_METADATA
- FASTQ_FILES_GUIDES
- WHITELIST 

## Outputs
-   guides_star_solo.h5ad


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
Your module workflow. You should import your processes here and structure the execution sequence



# Developer
```
docker run -it --platform linux/amd64 -v /Users/lf588/Downloads/pipeline_guide_counting/CRISPR-JAMBOREE/:/mnt/CRISPR-JAMBOREE ubuntu /bin/bash  
apt-get update  
apt-get install wget  
mkdir -p ~/miniconda3  
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh  
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3  
rm -rf ~/miniconda3/miniconda.sh  

source ~/.bashrc 

conda update -n base conda -y 
conda install -n base conda-libmamba-solver -y
conda config --set solver libmamba #this make conda very fast
#install nextflow
conda install bioconda::nextflow


```



