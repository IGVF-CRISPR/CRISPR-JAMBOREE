# `sceptre` guide assignment Nextflow process

The file `processes/guide_assignment_sceptre.nf` defines the Nextflow process `guide_assignment_sceptre`, which can be integrated into the IGVF pipeline.

## Docker image

The process uses the Docker image [ekatsevi/sceptre-igvf:v0.1](https://hub.docker.com/r/ekatsevi/sceptre-igvf/tags).

## Inputs

The only input is the path to the input `MuData` object, which must be passed into the process.

## Parameters

There are no parameters to this Nextflow process.

## Outputs

The process writes an output `MuData` file called `mudata_out.h5mu`. 

## Running the process

To run this process, the workflow in `main.nf` can be invoked via 
```
nextflow run main.nf --mudata_fp data/gasperini_guide_assignment_input.h5mu
```