#' Analysis of differential pathway activity scores
#'
#' @param x
#' @param treatment
#' @param control
#' @param group_by
#' @param metadata
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' data("pbmc_small")
#'
#' pbmc_small <- pbmc_small |>
#'   cat_score(gene_sets = "CP:KEGG", method = "GSVA")
#' pbmc_small
#'
#' diff_pathway <- compare_pathways(
#'   pbmc_small,
#'   control = "g2",
#'   treatment = "g1",
#'   group_by = "groups"
#' )
#' }
compare_pathways <- function(
    x,
    treatment,
    control,
    group_by = NULL,
    metadata = NULL) {
  if (inherits(x = x, what = "Seurat")) {
    if (!is.null(x = group_by)) {
      Idents(object = x) <- group_by
      x <- subset(x = x, idents = c(control, treatment))
    }
    # select which data to use
    assay <- "score"
    data_matrix <- x[[assay]]@counts
    metadata <- x@meta.data
  } else {
    data_matrix <- x
    if (is.null(metadata)) {
      stop("Please provide the metadata data.frame,
           requiring that the row name of the metadata is the sample name,
           and then the column name contains the grouping information of
           each sample.")
      # TODO: such as
    }
    metadata <- metadata
  }
  # DE
  group <-
    factor(metadata[, group_by], levels = c(treatment, control))
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  rownames(design) <- colnames(data_matrix)
  contrasts <- paste0(treatment, "-", control)
  compare <- limma::makeContrasts(
    contrasts = contrasts,
    levels = design
  )
  fit <- limma::lmFit(data_matrix, design)
  fit2 <- limma::contrasts.fit(fit, compare)
  fit3 <- limma::eBayes(fit2)
  diff <- limma::topTable(fit3, coef = 1, number = Inf)
  diff <- diff |>
    tibble::rownames_to_column("term") |>
    dplyr::mutate(compare = paste0(treatment, "_vs_", control))
  return(diff)
}
