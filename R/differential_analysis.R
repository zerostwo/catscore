#' Analysis of differential pathway activity scores
#'
#' @param x Pathway activity scores which can be given either as a Seurat
#' object, or as a matrix of expression values where rows correspond to
#' pathways and columns correspond to samples.
#' @param treatment Treatment group name.
#' @param control Control group name
#' @param group_by Column names in the metadata that contains grouping
#' information for the treatment and control groups.
#' @param metadata Metadata for each sample. Should be a data.frame where the
#' rows are sample names and the columns are additional metadata fields. Row
#' names in the metadata need to match the column names of the
#' pathway activity scores.
#'
#' @return A data.frame of differential pathway activity scores.
#'
#' @importFrom stats model.matrix
#'
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
      Seurat::Idents(object = x) <- group_by
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
  if (!package_check("limma", error = FALSE)) {
    stop("Please install limma:
         https://bioconductor.org/packages/release/bioc/html/limma.html")
  }
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
  database <- stringr::str_split_fixed(
    string = rownames(diff),
    pattern = "_", n = 2
  )[, 1]
  pathways <- stringr::str_split_fixed(
    string = rownames(diff),
    pattern = "_", n = 2
  )[, 2] |>
    stringr::str_replace_all(
      pattern = "_", replacement = " "
    ) |>
    stringr::str_to_sentence()
  diff <- diff |>
    tibble::rownames_to_column("term") |>
    dplyr::mutate(
      database = database,
      pathway = pathways,
      .after = "term"
    ) |>
    dplyr::mutate(
      compare = paste0(
        treatment, "_vs_", control
      )
    )
  return(diff)
}
