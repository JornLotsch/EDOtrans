#Function to create a probability matrix based on multimodal Gaussian distribution
#' @importFrom methods hasArg
#' @importFrom stats dnorm
CreateGMM <- function(Means, SDs, Weights, n = 1000, Prob = FALSE) {
  if (!hasArg("Means") | !hasArg("SDs") | !hasArg("Weights"))
    stop("CreateGMM: Incomplete parameters.")
  if (length(c(Means, SDs, Weights)) %% 3 != 0)
    stop("CreateGMM: Unequal number of modes in parameters.")
  sumWeights <- sum(Weights)
  if (sumWeights != 1) {
    Weights <- Weights / sumWeights
    warning("GMMInnerInterDistances: Weigthts changed to sum up to 1.",
            call. = FALSE)
  }

  if (Prob == FALSE) {
    GMMparam <- rbind(Weights * n, Means, SDs)
    GMMparamRO <- split(GMMparam, rep(1:ncol(GMMparam), each = nrow(GMMparam)))
    DataDF <- cbind.data.frame(Data = unlist(lapply(GMMparamRO, function(x) {
      do.call(rnorm, as.list(x))
    })),
      Cls = rep(1:length(Weights), unlist(lapply(GMMparamRO, "[[", 1))))
  } else {
    rangeX <- range(c(Means - 2 * SDs, Means + 2 * SDs))
    x <- seq(from = rangeX[1], to = rangeX[2], length.out = n)
    DataDF_wide <- data.frame(mapply(
        function(w, mean, sd)
          w * dnorm(x, mean, sd),
        mean = Means,
        sd = SDs,
        w = Weights
       ))
    DataDF <- cbind.data.frame(
      Data = rep(x, length(Means)),
      Prob = as.vector(as.matrix(DataDF_wide)),
      Cls = rep(c(1:length(Means)), each = length(x))
    )
  }

  #Return results
  return(DataDF)
}
