


#this is to speedup the conda solver
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
```

```bash

#Edit the input.config config file to match the directory that you extracted your data.



```
```bash
#Execute the nextflow from the project directory

nextflow main.nf -config input.config -resume


```