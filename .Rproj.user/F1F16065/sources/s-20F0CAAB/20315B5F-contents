#Cuts a vector into separate bfoups according to input breacks
#' @importFrom pracma size

cutGMM <- function(x, breaks, right = TRUE) {
  if (pracma::size(x)[2] > 1 | length(x) == 0)
    stop("CutGMM: x should be a single vector of length > 0.")
  if (length(breaks) > 0) {
    if (hasArg(right) == FALSE) {
      right = TRUE
    }
    breaks = sort(breaks)
    if(right == TRUE) {
      if(length(breaks) > 1) {
        compMat <- vapply(x, function(x) x > breaks, logical(length(breaks)))
        ClassesGMM <- colSums(compMat)+1
      } else ClassesGMM <- as.numeric(x > breaks) + 1
    } else {
      if(length(breaks) > 1) {
        compMat <- vapply(x, function(x) x >= breaks, logical(length(breaks)))
        ClassesGMM <- colSums(compMat)+1
      } else ClassesGMM <- as.numeric(x >= breaks) + 1
    }
  } else {
    warning("cutGMM: No breaks provided. Assuming one single class.", call. = FALSE)
  }
  return(ClassesGMM)
}
