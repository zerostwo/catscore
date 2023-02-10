#' Construct gene sets
#'
#' @param gene_sets MSigDB collection abbreviation, such as H or CP:KEGG.
#' @param species Species name, such as human or mouse.
#'
#' @return A list of vectors of genes for biological pathway.
#'
#' @importFrom msigdbr msigdbr_collections msigdbr
#'
#' @export
#'
#' @examples
#' gene_sest <- construct_gene_sets(gene_sets = "H", species = "human")
construct_gene_sets <- function(gene_sets,
                                species = "human") {
  # TODO: Gene sets from GO, KEGG, and RECTOME
  # Gene sets from msigdbr
  collections <- msigdbr_collections()
  gs_cat <- names(table(collections$gs_cat))
  gs_subcat <- names(table(collections$gs_subcat))[-1]
  if (all(gene_sets %in% gs_cat)) {
    gene_sets <- lapply(gene_sets, function(x) {
      msigdbr(
        species = species,
        category = x
      )[, c("gs_name", "gene_symbol")]
    })
    gene_sets <- do.call(rbind, gene_sets)
  } else if (all(gene_sets %in% gs_subcat)) {
    gene_sets <- lapply(gene_sets, function(x) {
      msigdbr(
        species = species,
        subcategory = x
      )[, c("gs_name", "gene_symbol")]
    })
    gene_sets <- do.call(rbind, gene_sets)
  } else {
    stop("The gene sets name is wrong!")
  }
  gene_sets_list <- split(gene_sets$gene_symbol, gene_sets$gs_name)
  return(gene_sets_list)
}

#' Pathway activity scores (PASs)
#'
#' @param x Gene expression data which can be given either as a Seurat object,
#' or as a matrix of expression values where rows correspond to genes and
#' columns correspond to samples. This matrix can be also in a sparse format,
#' as a dgCMatrix.
#' @param gene_sets A list of vectors of genes for biological pathway.
#' @param species Species name, such as human or mouse.
#' @param method Method to employ in estimation of pathway activity scores
#' (PASs). By default this is set to AUCell (Sara Aibar et al., Nature
#' Methods, 2017) and other options are VISION (David DeTomaso et al.,
#' Nature Communications, 2019), ssGSEA (David A. Barbie et al., Nature, 2009),
#' or GSVA (Sonja Hänzelmann et al., BMC Bioinformatics, 2013).
#' @param ncores Number of cores to use for computation.
#' @param return_df Whether to return a data.frame.
#' @param ... Other arguments.
#'
#' @return A data.frame of pathway activity scores.
#'
#' @importFrom Matrix t as.matrix
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("pbmc_small")
#' score_data <- pbmc_small |>
#'   cat_score(
#'     gene_sets = "CP:KEGG",
#'     method = "gsva"
#'   )
#' }
cat_score <- function(x,
                      gene_sets,
                      species = "human",
                      method = "AUCell",
                      ncores = 1,
                      return_df = FALSE,
                      ...) {
  match.arg(method, c("AUCell", "GSVA", "ssGSEA", "VISION"))
  # Get counts
  if (inherits(x = x, what = "Seurat")) {
    assay <- "RNA"
    slot <- "counts"
    data_matrix <-
      as.matrix(Seurat::GetAssayData(x, assay = assay, slot = slot))
  } else {
    data_matrix <- as.matrix(x)
  }
  data_matrix <- data_matrix[rowSums(data_matrix) > 0, ]

  # Construct gene sets
  gene_sets <- construct_gene_sets(
    gene_sets = gene_sets,
    species = species
  )

  # Score
  # VISION
  if (method == "VISION") {
    if (!package_check("VISION", error = FALSE)) {
      stop("Please install VISION: https://github.com/YosefLab/VISION")
    }
    signatures <- purrr::map(names(gene_sets), function(x) {
      sig_data <- rep(1, length(gene_sets[[x]]))
      names(sig_data) <- gene_sets[[x]]
      VISION::createGeneSignature(
        name = x,
        sigData = sig_data
      )
    })
    score_res <- VISION::Vision(
      data = methods::as(data_matrix, "dMatrix"),
      signatures = signatures,
      ...
    )
    options(mc.cores = ncores)
    score_res <- VISION::analyze(score_res)
    score_data <- as.data.frame(t(score_res@SigScores))
    colnames(score_data) <- colnames(data_matrix)
  }
  # AUCelll
  if (method == "AUCell") {
    if (!package_check("AUCell", error = FALSE)) {
      stop("Please install AUCell:
           https://bioconductor.org/packages/release/bioc/html/AUCell.html")
    }
    # 1. Build gene-expression rankings for each cell
    cells_rankings <- AUCell::AUCell_buildRankings(
      exprMat = data_matrix,
      plotStats = FALSE,
      nCores = ncores,
      ...
    )

    # 2. Calculate enrichment for the gene signatures (AUC)
    score_raw <-
      AUCell::AUCell_calcAUC(
        geneSets = gene_sets,
        rankings = cells_rankings,
        nCores = ncores,
        ...
      )
    score_data <- as.data.frame(AUCell::getAUC(score_raw))
  }
  # GSVA
  if (method == "GSVA") {
    if (!package_check("GSVA", error = FALSE)) {
      stop("Please install GSVA:
           https://bioconductor.org/packages/release/bioc/html/GSVA.html")
    }
    score_data <-
      GSVA::gsva(
        expr = data_matrix,
        gset.idx.list = gene_sets,
        method = c("gsva"),
        kcdf = c("Poisson"),
        tau = 1,
        min.sz = 2,
        parallel.sz = ncores,
        BPPARAM = BiocParallel::SerialParam(progressbar = TRUE),
        ...
      )
    score_data <- as.data.frame(score_data)
  }
  # ssGSEA
  if (method == "ssGSEA") {
    if (!package_check("GSVA", error = FALSE)) {
      stop("Please install GSVA:
           https://bioconductor.org/packages/release/bioc/html/GSVA.html")
    }
    score_data <-
      GSVA::gsva(
        expr = data_matrix,
        gset.idx.list = gene_sets,
        method = c("ssgsea"),
        kcdf = c("Poisson"),
        tau = 0.25,
        min.sz = 2,
        ssgsea.norm = TRUE,
        parallel.sz = ncores,
        BPPARAM = BiocParallel::SerialParam(progressbar = TRUE),
        ...
      )
    score_data <- as.data.frame(score_data)
  }
  if (return_df) {
    return(score_data)
  } else if (inherits(x = x, what = "Seurat")) {
    x[["score"]] <- create_assay_object(counts = score_data)
    return(x)
  }
}

create_assay_object <-
  function(counts) {
    if (!missing(x = counts)) {
      if (anyDuplicated(x = rownames(x = counts))) {
        warning(
          "Non-unique features (rownames) present in the input matrix,
          making unique",
          call. = FALSE,
          immediate. = TRUE
        )
        rownames(x = counts) <-
          make.unique(names = rownames(x = counts))
      }
      if (anyDuplicated(x = colnames(x = counts))) {
        warning(
          "Non-unique cell names (colnames) present in the input matrix,
          making unique",
          call. = FALSE,
          immediate. = TRUE
        )
        colnames(x = counts) <-
          make.unique(names = colnames(x = counts))
      }
      if (is.null(x = colnames(x = counts))) {
        stop("No cell names (colnames) names present in the input matrix")
      }
      if (any(rownames(x = counts) == "")) {
        stop("Feature names of counts matrix cannot be empty",
          call. = FALSE
        )
      }
      if (nrow(x = counts) > 0 && is.null(x = rownames(x = counts))) {
        stop("No feature names (rownames) names present in the input matrix")
      }
      if (!inherits(x = counts, what = "dgCMatrix")) {
        if (inherits(x = counts, what = "data.frame")) {
          counts <- Seurat::as.sparse(x = counts)
        }
      }
      data <- counts
    }
    if (!is.vector(x = rownames(x = counts))) {
      rownames(x = counts) <- as.vector(x = rownames(x = counts))
    }
    if (!is.vector(x = colnames(x = counts))) {
      colnames(x = counts) <- as.vector(x = colnames(x = counts))
    }
    if (!is.vector(x = rownames(x = data))) {
      rownames(x = data) <- as.vector(x = rownames(x = data))
    }
    if (!is.vector(x = colnames(x = data))) {
      colnames(x = data) <- as.vector(x = colnames(x = data))
    }
    init_meta_features <- data.frame(row.names = rownames(x = data))
    assay <- methods::new(
      Class = "Assay",
      counts = counts,
      data = data,
      scale.data = methods::new(Class = "matrix"),
      meta.features = init_meta_features,
      misc = list()
    )
    return(assay)
  }
