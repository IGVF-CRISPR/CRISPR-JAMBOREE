# `sceptre` inference Nextflow process

The file `main.nf` defines the Nextflow process `inference_sceptre`, which can be integrated into the IGVF pipeline.

## Docker image

The process uses the Docker image [ekatsevi/sceptre-igvf:v0.1](https://hub.docker.com/r/ekatsevi/sceptre-igvf/tags).

## Inputs

The only input is the path to the input `MuData` object, which must be passed into the process.

## Parameters

The parameters used by the process, which must be specified in `input.config`, are as follows:

- `side`
- `grna_integration_strategy`
- `resampling_approximation`
- `control_group`
- `resampling_mechanism`
- `formula_object`

These parameters are documented [here](https://timothy-barry.github.io/sceptre-book/set-analysis-parameters.html). Note that the `formula_object` should be passed as a string, such as `"~ log(response_n_nonzero) + log(response_n_umis)"`. All parameters have sensible default values.

## Outputs

The process writes an output `MuData` file called `mudata_out.h5mu`. 

## Running the process

To run this process, the workflow in `main.nf` can be invoked via 
```
nextflow run main.nf \
    --mudata_fp data/gasperini_inference_input.h5mu \
    -config input.config
```