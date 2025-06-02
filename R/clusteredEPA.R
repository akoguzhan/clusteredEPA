#' General Wrapper for Clustered EPA Tests
#'
#' @description
#' Runs one of several clustered EPA tests on a balanced panel data frame.
#'
#' @param df Balanced panel data.frame.
#' @param DL Character; name of the loss differential variable in df.
#' @param H Character vector; names of conditioning variables in df (optional).
#' @param id Character; name of the unit identifier column.
#' @param time Character; name of the time identifier column.
#' @param test Character; which test to run: "overall_EPA_test", "epa_clustered_known", "epa_clustered_split", or "epa_clustered_selective" (default).
#' @param ... Additional arguments passed to the selected test function.
#'
#' @return Output of the selected test function.
#' @export
clusteredEPA <- function(df, DL, H = NULL, id, time, test = "epa_clustered_selective", ...) {
  # Input checks
  if (!is.data.frame(df)) stop("df must be a data frame.")
  if (!(DL %in% names(df))) stop("Loss differential variable not found in df.")
  if (!(id %in% names(df))) stop("id column not found in df.")
  if (!(time %in% names(df))) stop("time column not found in df.")
  if (!is.null(H) && !all(H %in% names(df))) stop("Some conditioning variables not found in df.")

  # Build Z matrix: always include constant (loss differential), then DL*H for each H
  if (is.null(H)) {
    Z <- as.matrix(df[[DL]])
    colnames(Z) <- DL
  } else {
    Z <- sapply(H, function(h) df[[DL]] * df[[h]])
    Z <- cbind(df[[DL]], Z)
    colnames(Z) <- c(DL, paste0(DL, "_x_", H))
  }

  # Prepare id and time vectors
  id_vec <- df[[id]]
  time_vec <- df[[time]]

  # Dispatch to the selected test
  test <- match.arg(test, c("overall_EPA_test", "epa_clustered_known", "epa_clustered_split", "epa_clustered_selective"))

  if (test == "overall_EPA_test") {
    res <- overall_EPA_test(Z = Z, id = id_vec, time = time_vec, ...)
  } else if (test == "epa_clustered_known") {
    res <- epa_clustered_known(Z = Z, id = id_vec, time = time_vec, ...)
  } else if (test == "epa_clustered_split") {
    res <- epa_clustered_split(df = df, id = id, time = time, Z_names = colnames(Z), ...)
  } else if (test == "epa_clustered_selective") {
    res <- epa_clustered_selective(Z = Z, id = id_vec, time = time_vec, ...)
  } else {
    stop("Unknown test type.")
  }

  return(res)
}