
## Runing our starSolo nextflow module to count guides
This workflow will receive the whitelist, fastq and metadata and generate and anndata (cells x guides) using starSolo.s

```bash
#iNSTALL CONDA
conda update -n base conda -y 
conda install -n base conda-libmamba-solver -y
conda config --set solver libmamba #this make conda very fast
#install nextflow
conda install bioconda::nextflowå

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