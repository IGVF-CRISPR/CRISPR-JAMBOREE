## Estimate dispersions for all genes in muData object using DESeq2's estimateDispersions function.

# required packages
library(MuData)
library(SummarizedExperiment)
library(DESeq2)


## Define functions --------------------------------------------------------------------------------

#' Fit negative binomial distributions using DESeq2
#' 
#' Fit negative binomial distributions to estimate dispersion for every gene in a
#' muData object.
#' 
#' @param mu A muData object containing scRNA-seq counts as SummarizedExperiment.
#' @param size_factors Different options to compute size factors. 'libsize' are simple size factors
#'   based on library size, for others see DESeq2 manual for more information.
#' @param fit_type Options for dispersion fit method, see DESeq2 manual for more information.
#' @param disp_type Type of dispersion to return and add to muData object, see DESeq2 manual for 
#'   more information.
#'   
#' @return A muData object with estimated dispersions per gene in rowData of scRNA-seq slot.
fit_negbinom_deseq2 <- function(mu,
                                size_factors = c("ratio", "poscounts", "iterate", "libsize"),
                                fit_type = c("parametric", "local", "mean"),
                                disp_type = c("dispersion", "dispFit", "dispGeneEst", "dispMAP")) {
  
  # process optional arguments and get default values
  size_factors <- match.arg(size_factors)
  fit_type <- match.arg(fit_type)
  disp_type <- match.arg(disp_type)
  
  # get SummarizedExperiment containing scRNA-seq data
  rna_se <- mu@ExperimentList$gene
  
  # check if RNA SE object already contains dispersion and raise warning if data will be overwritten
  present_row_data <- colnames(rowData(rna_se)) %in% c("mean", "dispersion", "disp_outlier_deseq2")
  present_col_data <- colnames(colData(rna_se)) == "size_factors"
  if (any(present_row_data, present_col_data)) {
    warning("Dispersion data found in muData object, will overwrite values", call. = FALSE)
  }
  
  # create DESeq2 object containing count data
  dds <- DESeqDataSetFromMatrix(countData = assay(rna_se, 1),
                                colData = colData(rna_se),
                                design = ~ 1)
  
  # compute size factors
  if (size_factors == "libsize") {
    total_umis <- colSums(assay(dds, assay))
    manual_size_factors <- total_umis / mean(total_umis)
    sizeFactors(dds) <- manual_size_factors
  } else {
    dds <- estimateSizeFactors(dds, type = size_factors)
  }
  
  # estimate dispersion
  dds <- estimateDispersions(dds, fitType = fit_type)
  
  # add mean expression, dispersion and cell size factors to rowData and colData of RNA SE object
  rowData(rna_se)[, "mean"] <- rowData(dds)[, "baseMean"]
  rowData(rna_se)[, "dispersion"] <- rowData(dds)[, disp_type]
  rowData(rna_se)[, "disp_outlier_deseq2"] <- rowData(dds)[, "dispOutlier"]
  colData(rna_se)[, "size_factors"] <- sizeFactors(dds)
  
  # store dispersion function in metadata of RNA SE object
  metadata(rna_se)[["dispersionFunction"]] <- dispersionFunction(dds)
  
  # add RNA SE with estimated dispersion back into muData object
  mu@ExperimentList$gene <- rna_se
  
  return(mu)
  
}


## Estimate dispersions for input muData object ----------------------------------------------------

## Replace with NextFlow input and output
infile <- "output/gasperini_simulation_input.h5mu"
outfile <- paste0(tools::file_path_sans_ext(infile), "_disp.h5mu")

# load muData object
mu <- readH5MU(infile)

# estimate dispersion and add estimates to muData object
mu <- fit_negbinom_deseq2(mu, size_factors = "poscounts", fit_type = "parametric",
                          disp_type = "dispersion")

# write to output file
writeH5MU(mu, file = outfile)
