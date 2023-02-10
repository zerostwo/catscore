#' Construct gene sets
#'
#' @param gene_sets
#' @param species
#'
#' @return
#' @export
#'
#' @examples
construct_gene_sets <- function(gene_sets,
                                species = "human") {
  # TODO: Gene sets from GO, KEGG, and RECTOME
  # Gene sets from msigdbr
  collections <- msigdbr::msigdbr_collections()
  gs_cat <- names(table(collections$gs_cat))
  gs_subcat <- names(table(collections$gs_subcat))[-1]
  if (all(gene_sets %in% gs_cat)) {
    gene_sets <- purrr::map_dfr(gene_sets, function(x) {
      msigdbr::msigdbr(
        species = species,
        category = x
      )[, c("gs_name", "gene_symbol")]
    })
  } else if (all(gene_sets %in% gs_subcat)) {
    gene_sets <- purrr::map_dfr(gene_sets, function(x) {
      msigdbr::msigdbr(
        species = species,
        subcategory = x
      )[, c("gs_name", "gene_symbol")]
    })
  } else {
    stop("The gene sets name is wrong!")
  }
  gene_sets_list <- split(gene_sets$gene_symbol, gene_sets$gs_name)
  return(gene_sets_list)
}

#' Pathway activity scores (PASs)
#'
#' @param x
#' @param gene_sets
#' @param species
#' @param method
#' @param ncores
#' @param ...
#'
#' @return
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
                      return_seurat_obj = TRUE,
                      ...) {
  match.arg(method, c("AUCell", "GSVA", "ssGSEA", "VISION"))
  # Get counts
  if (inherits(x = x, what = "Seurat")) {
    assay <- "RNA"
    slot <- "counts"
    data_matrix <-
      Matrix::as.matrix(Seurat::GetAssayData(x, assay = assay, slot = slot))
  } else {
    data_matrix <- Matrix::as.matrix(x)
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
    signatures <- purrr::map(names(gene_sets), function(x) {
      sig_data <- rep(1, length(gene_sets[[x]]))
      names(sig_data) <- gene_sets[[x]]
      VISION::createGeneSignature(
        name = x,
        sigData = sig_data
      )
    })
    score_res <- VISION::Vision(
      data = as(data_matrix, "dMatrix"),
      signatures = signatures,
      ...
    )
    options(mc.cores = ncores)
    score_res <- VISION::analyze(score_res)
    score_data <- as.data.frame(Matrix::t(score_res@SigScores))
    colnames(score_data) <- colnames(data_matrix)
  }
  # AUCelll
  if (method == "AUCell") {
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

  if (inherits(x = x, what = "Seurat") && return_seurat_obj) {
    x[["score"]] <- create_assay_object(counts = score_data)
    return(x)
  } else {
    return(score_data)
  }
}

create_assay_object <-
  function(counts,
           data,
           check_matrix = FALSE) {
    if (missing(x = counts) && missing(x = data)) {
      stop("Must provide either 'counts' or 'data'")
    } else if (!missing(x = counts) && !missing(x = data)) {
      stop("Either 'counts' or 'data' must be missing; both cannot be provided")
    } else if (!missing(x = counts)) {
      if (anyDuplicated(x = rownames(x = counts))) {
        warning(
          "Non-unique features (rownames) present in the input matrix, making unique",
          call. = FALSE,
          immediate. = TRUE
        )
        rownames(x = counts) <-
          make.unique(names = rownames(x = counts))
      }
      if (anyDuplicated(x = colnames(x = counts))) {
        warning(
          "Non-unique cell names (colnames) present in the input matrix, making unique",
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
      if (isTRUE(x = check_matrix)) {
        SeuratObject::CheckMatrix(object = counts)
      }
      data <- counts
    } else if (!missing(x = data)) {
      if (anyDuplicated(x = rownames(x = data))) {
        warning(
          "Non-unique features (rownames) present in the input matrix, making unique",
          call. = FALSE,
          immediate. = TRUE
        )
        rownames(x = data) <- make.unique(names = rownames(x = data))
      }
      if (anyDuplicated(x = colnames(x = data))) {
        warning(
          "Non-unique cell names (colnames) present in the input matrix, making unique",
          call. = FALSE,
          immediate. = TRUE
        )
        colnames(x = data) <- make.unique(names = colnames(x = data))
      }
      if (is.null(x = colnames(x = data))) {
        stop("No cell names (colnames) names present in the input matrix")
      }
      if (any(rownames(x = data) == "")) {
        stop("Feature names of data matrix cannot be empty",
          call. = FALSE
        )
      }
      if (nrow(x = data) > 0 && is.null(x = rownames(x = data))) {
        stop("No feature names (rownames) names present in the input matrix")
      }
      counts <- new(Class = "matrix")
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
    assay <- new(
      Class = "Assay",
      counts = counts,
      data = data,
      scale.data = new(Class = "matrix"),
      meta.features = init_meta_features,
      misc = list()
    )
    return(assay)
  }
