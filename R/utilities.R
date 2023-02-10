package_check <- function(..., error = TRUE) {
  pkgs <- unlist(x = c(...), use.names = FALSE)
  package_installed <- vapply(
    X = pkgs,
    FUN = requireNamespace,
    FUN.VALUE = logical(length = 1L),
    quietly = TRUE
  )
  if (error && any(!package_installed)) {
    stop(
      "Cannot find the following packages: ",
      paste(pkgs[!package_installed], collapse = ", "),
      ". Please install"
    )
  }
  invisible(x = package_installed)
}
