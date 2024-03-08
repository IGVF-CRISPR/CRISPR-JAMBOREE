#' Negative Binomial Test
#'
#' @param mudata_input_fp Path to input MuData
#' @param mudata_output_fp Path to output MuData
#'
#' @return This function is called for its side effect of writing the output MuData.
# TODO: confirm that Low MOI and high MOI case are the same for negative binomial
perform_negbinom_regression <- function(mudata_input_fp, mudata_output_fp) {
    
  # read in input mudata
  mudata <- MuData::readH5MU(mudata_input_fp)

  # Rename primary assay to 'counts' for both guide and gene
  if(is.null(SummarizedExperiment::assayNames(mudata[["gene"]]))){
    SummarizedExperiment::assayNames(mudata[['gene']]) <- 'counts'    
  } else{
    SummarizedExperiment::assayNames(mudata[['gene']])[[1]] <- 'counts'  
  }
  SummarizedExperiment::assayNames(mudata[['guide']])[[1]] <- 'counts'

  # get MOI
  moi <- MultiAssayExperiment::metadata(mudata[["guide"]])$moi

  # In low-MOI case, extract control cells as those containing an NT gRNA
  if (moi == "low") {
    non_targeting_guides <- SummarizedExperiment::rowData(mudata[["guide"]]) |>
      as.data.frame() |>
      dplyr::filter(targeting == "FALSE") |>
      rownames()
    nt_grna_presence <- SummarizedExperiment::assay(
      mudata[["guide"]],
      "guide_assignment"
    )[non_targeting_guides, ] |>
      apply(MARGIN = 2, FUN = max)
    control_cells <- names(nt_grna_presence)[nt_grna_presence == 1]
  }

  # Extract pairs to test
  pairs_to_test <- MultiAssayExperiment::metadata(mudata)$pairs_to_test |> 
    as.data.frame()

  # Initialize test results data frame based on pairs_to_test
  test_results <- pairs_to_test |>
    dplyr::mutate(p_value = NA_real_, l2fc = NA_real_)
    
  # Carry out the negative binomial regression test for each pair
  for (pair_idx in 1:nrow(pairs_to_test)) {
        
    # Extract the gene and element to be tested
    gene_id <- pairs_to_test[pair_idx, "gene_id"]
    intended_target_name <- pairs_to_test[pair_idx, "intended_target_name"]

    # find vector with whether element is targeted or not
    grnas_targeting_element <- SummarizedExperiment::rowData(mudata[["guide"]]) |>
      as.data.frame() |>
    dplyr::filter(intended_target_name == !!intended_target_name) |>
      rownames()
    element_targeted <- SummarizedExperiment::assay(
      mudata[["guide"]],
      "guide_assignment"
    )[grnas_targeting_element, ] |>
      apply(MARGIN = 2, FUN = max)
    
    # define treatment cells for low MOI case
    treatment_cells <- names(element_targeted)[element_targeted == 1]

    # get gene expression for gene of interest
    gene_expression = SummarizedExperiment::assay(
      mudata[["gene"]],
      "counts"
    )[gene_id, ]

    # get total gene expression for each cell
    umis_per_cell <- SummarizedExperiment::colData(mudata[["gene"]])
    umis_per_cell <- umis_per_cell[, 'total_gene_umis']
        
    # create modeling data frame
    model_df <- data.frame(cbind(element_targeted, gene_expression, umis_per_cell))

    # subset cells in low MOI case - cells with guide and NT guide
    if (moi == "low") {
      model_df <- model_df[c(control_cells, treatment_cells), ]
    }

    print(dim(model_df))
    
    # run negative binomial regression
    mdl <- MASS::glm.nb(gene_expression ~ element_targeted + offset(log(umis_per_cell)), model_df)

    # save L2FC and p-value
    mdl.coeffs <- summary(mdl)$coefficients
    test_results[pair_idx, "p_value"] <- mdl.coeffs['element_targeted', 'Pr(>|z|)']
    mdl.intercept <- mdl.coeffs['(Intercept)', 'Estimate']
    mdl.beta1 <- mdl.coeffs['element_targeted', 'Estimate']
    mdl.l2fc <- exp(mdl.intercept + mdl.beta1) / exp(mdl.intercept)
    test_results[pair_idx, 'l2fc'] <- mdl.l2fc
  }

  # Add output to MuData and write to disk
  mudata_output <- mudata
  MultiAssayExperiment::metadata(mudata_output)$test_results <- test_results
  MuData::writeH5MU(mudata_output, mudata_output_fp)
}
       