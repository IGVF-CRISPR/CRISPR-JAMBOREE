# Power analysis/simulation

## Install R dependencies
To run the simulation engine on the jamboree cluster, install required R dependencies in a personal library.

Create a personal library in your home directory using the following commands in a terminal:
```
# create personal linrary directory
mkdir "~/R"

# create .Renviron file pointing to personal library
echo "R_LIBS_USER=~/R" > ~/.Renviron
```

Launch an R console and run this command to check that your personal library is used by R:
```
.libPaths()
```

Now you can install the required packages using this command:
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
```
