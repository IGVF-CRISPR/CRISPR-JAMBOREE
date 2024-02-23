#' Wilcoxon test
#'
#' @param mudata_input_fp Path to input MuData
#' @param mudata_output_fp Path to output MuData
#' @param side The sidedness of the test (`left`, `right`, or `both`)
#'
#' @return This function is called for its side effect of writing the output MuData.
compute_wilcoxon_test <- function(mudata_input_fp, mudata_output_fp, side) {
  # Read input MuData
  mudata <- MuData::readH5MU(mudata_input_fp)
  # Rename primary assay to 'counts'
  if(is.null(SummarizedExperiment::assayNames(mudata[['gene']]))){
    SummarizedExperiment::assayNames(mudata[['gene']]) <- 'counts'    
  } else{
    SummarizedExperiment::assayNames(mudata[['gene']])[[1]] <- 'counts'  
  }
  SummarizedExperiment::assayNames(mudata[['guide']])[[1]] <- 'counts'
  # Extract pairs to test and MOI
  pairs_to_test <- MultiAssayExperiment::metadata(mudata)$pairs_to_test |> 
    as.data.frame()
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
  # Initialize test results data frame based on pairs_to_test
  test_results <- pairs_to_test |>
    dplyr::mutate(p_value = NA_real_)
  # Carry out the Wilcoxon test for each pair
  for (pair_idx in 1:nrow(pairs_to_test)) {
    # Extract the gene and element to be tested
    gene_id <- pairs_to_test[pair_idx, "gene_id"]
    intended_target_name <- pairs_to_test[pair_idx, "intended_target_name"]
    # Find the set of treatment cells (i.e. cells with element targeted)
    grnas_targeting_element <- SummarizedExperiment::rowData(mudata[["guide"]]) |>
      as.data.frame() |>
      dplyr::filter(intended_target_name == !!intended_target_name) |>
      rownames()
    element_targeted <- assay(
      mudata[["guide"]],
      "guide_assignment"
    )[grnas_targeting_element, ] |>
      apply(MARGIN = 2, FUN = max)
    treatment_cells <- names(element_targeted)[element_targeted == 1]
    # Set controls cells using the complement set in high MOI
    if (moi == "high") {
      control_cells <- names(element_targeted)[element_targeted == 0]
    }
    # extract expressions for treatment and control cells
    treatment_expressions <- SummarizedExperiment::assay(
      mudata[["gene"]],
      "counts"
    )[gene_id, treatment_cells]
    control_expressions <- SummarizedExperiment::assay(
      mudata[["gene"]],
      "counts"
    )[gene_id, control_cells]
    # Map `side` argument to `alternative` argument required by `wilcox.test()`
    alternative <- switch(side,
      left = "less",
      right = "greater",
      both = "two.sided"
    )
    # Carry out the Wilcoxon test
    test_results[pair_idx, "p_value"] <- stats::wilcox.test(
      x = treatment_expressions,
      y = control_expressions,
      alternative = alternative
    )$p.value
  }
  # Add output to MuData and write to disk
  mudata_output <- mudata
  MultiAssayExperiment::metadata(mudata_output)$test_results <- test_results
  MuData::writeH5MU(mudata_output, mudata_output_fp)
}