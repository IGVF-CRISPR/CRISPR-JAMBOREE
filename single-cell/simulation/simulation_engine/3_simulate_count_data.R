## Simulate Perturb-seq data from a muData object containing perturbation effects to simulate and
## dispersion estimates per gene

# required packages
library(MuData)
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(tibble)


## Define functions --------------------------------------------------------------------------------

#' Simulate Perturb-seq data
#' 
#' Simulate Perturb-seq data from a muData object containing perturbation effects to simulate and
#' dispersion estimates per gene.
#' 
#' @param mu A muData object containing scRNA-seq counts, perturbation status and simulation
#'   parameters.
#' @param guide_var Degree of guide-to-guide variability in effect sizes. Used to draw effect sizes
#'   for individual guides from a normal distribution with mean = specified effect size and
#'   sd = guide_var.
#'   
#' @return A muData object containing simulated Perturb-seq data.
simulate_pert_data <- function(mu, guide_var = 0) {
  
  # get pairs to test table containing pairs to perturb
  pairs_to_test <- data.frame(mu@metadata$pairs_to_test, stringsAsFactors = FALSE)
  
  # get perturbation status data from muData object
  pert_status <- assay(mu@ExperimentList$guide, "guide_assignment")
  pert_meta <- rowData(mu@ExperimentList$guide)
  
  # create gene x cell effect size matrix
  effect_size_matrix <- create_effect_size_mat(mu@ExperimentList$gene, pairs_to_test, pert_status,
                                               pert_meta, guide_var)
  
  # simulate scRNA-seq counts using dispersion and mean stored in muData object with perturbation
  # effect sizes in the effect size matrix
  sim_counts <- simulate_counts(mu@ExperimentList$gene, effect_size_matrix)
  
  # create SummarizedExperiment object for simulated scRNA-seq counts
  sim_gene <- mu@ExperimentList$gene
  assay(sim_gene) <- sim_counts
  
  # add to muData object
  mu@ExperimentList$gene <- sim_gene

  return(mu)
  
}

# function to create perturbation effects per perturbed cell
create_pert_effect_sizes <- function(pairs_to_test, pert_status, pert_meta, guide_var) {
  
  # create table with guides per perturbation - gene pair to test
  pert_meta$guide <- rownames(pert_meta)
  guide_targets <- pert_meta[, c("guide", "intended_target_name")]
  pairs_to_test <- merge(pairs_to_test, guide_targets, by = "intended_target_name")
  
  # pick guide level effect sizes based on target level effect size and guide variability
  guide_var <- ifelse(pairs_to_test$effect_size != 1, guide_var, 0)
  pairs_to_test$guide_effect_size <- mapply(pairs_to_test$effect_size, guide_var,
                                            FUN = rnorm, MoreArgs = list(n = 1))
  
  # get all guides per cell and add to pairs_to_test to create a table of guide effects per cell
  # in long format
  guides_per_cell <- get_guides_per_cell(pert_status)
  guide_effects_per_cell <- merge(guides_per_cell, pairs_to_test, by = "guide")
  
  # if a cell has more than one guide for a given target - gene pair, randomly pick one
  guide_effects_per_cell <- guide_effects_per_cell %>% 
    as.data.frame() %>% 
    group_by(cell, gene_id) %>% 
    slice_sample(n = 1)
  
  # select columns for output
  guide_effects_per_cell <- select(guide_effects_per_cell, cell, gene_id, guide_effect_size)
  
  return(guide_effects_per_cell)
  
}

# function to create gene-by-cells effect size matrix
create_effect_size_mat <- function(genex, pairs_to_test, pert_status, pert_meta, guide_var) {
  
  # create effect sizes for all perturbed cell-gene pairs
  effects_per_cell <- create_pert_effect_sizes(pairs_to_test, pert_status, pert_meta, guide_var)
  
  # create table with all cell-gene combinations
  cells_by_genes <- expand.grid(cell = colnames(genex), gene_id = rownames(genex))
  
  # add effect sizes for perturbed cell-gene pairs
  cells_by_genes <- merge(cells_by_genes, effects_per_cell, by = c("cell", "gene_id"), all.x = TRUE)
  cells_by_genes[is.na(cells_by_genes$guide_effect_size), "guide_effect_size"] <- 1
  
  # convert to gene-by-cells effect size matrix
  effect_size_matrix <- cells_by_genes %>% 
    pivot_wider(names_from = cell, values_from = guide_effect_size, values_fill = 1) %>% 
    column_to_rownames(var = "gene_id") %>% 
    as.matrix()
  
  return(effect_size_matrix)
  
}

# get a table with all guides detected in each cell
get_guides_per_cell <- function(pert_status) {
  
  # get list with guides per cell
  guides_per_cell <- apply(pert_status, MARGIN = 2, FUN = function(x) { which(x != 0) } )
  
  # get unique cell barcodes and repeat each cell barcode to match the number of guides per cell
  cells <- vapply(guides_per_cell, FUN = length, FUN.VALUE = integer(1))
  cells <- rep(names(cells), times = cells)
  
  # get vector of guides
  guides <- rownames(pert_status)[unlist(guides_per_cell)]
  
  # create table with all guides per cell barcode
  guides_per_cell_table <- data.frame(cell = cells, guide = guides)
  
  return(guides_per_cell_table)
  
}

# simulate scRNA-seq counts based on provided effect size matrix
simulate_counts <- function(genex, effect_size_matrix) {
  
  # get mean expression and dispersions for each gene in effect size matrix
  gene_means <- structure(rowData(genex)[["mean"]], names = rownames(genex))
  gene_dispersions <- structure(rowData(genex)[["dispersion"]], names = rownames(genex))
  gene_means <- gene_means[rownames(effect_size_matrix)]
  gene_dispersions <- gene_dispersions[rownames(effect_size_matrix)]
  
  # get size factors for all cells in effect size matrix
  cell_size_factors <- structure(colData(genex)[["size_factors"]], names = colnames(genex))
  cell_size_factors <- cell_size_factors[colnames(effect_size_matrix)]
  
  # number of genes and cells
  n_genes <- length(gene_means)
  n_cells <- length(cell_size_factors)
  
  # make mu matrix for simulation
  mu <- matrix(rep(gene_means, n_cells), ncol = n_cells)
  mu <- sweep(mu, 2, cell_size_factors, "*")  # add cell-to-cell variability based on size factors
  
  # inject perturbation effects by element-wise product of mu and effect_size_mat
  mu <- mu * effect_size_matrix
  
  # simulate counts
  sim_counts <- Matrix(rnbinom(n_cells * n_genes, mu = mu, size = 1 / gene_dispersions),
                       ncol = n_cells, dimnames = list(names(gene_means), names(cell_size_factors)))
  
  return(sim_counts)
  
}

## Run simulation ----------------------------------------------------------------------------------

## Replace with NextFlow input and output
infile <- "output/gasperini_simulation_input_disp.h5mu"
outfile <- "output/gasperini_simulation_output.h5mu"

# load muData file
mu <- readH5MU(infile)

# simulate Perturb-seq data
output <- simulate_pert_data(mu, guide_var = 0.1)

# write to output file
writeH5MU(mu, file = outfile)
