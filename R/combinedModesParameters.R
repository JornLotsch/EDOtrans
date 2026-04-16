# Calculates parameters of combined Gaussian modes
combinedModesParameters <- function(Means, SDs, Weights, n) {
  # Validation
  if (length(Means) != length(SDs) || length(Means) != length(Weights)) {
    stop("Means, SDs, and Weights must have the same length")
  }
  if (length(n) != 1 || !is.numeric(n) || n <= 0) {
    stop("n must be a positive number")
  }

  # Order by means
  Means0 <- Means[order(Means)]
  SDs0 <- SDs[order(Means)]
  Weights0 <- Weights[order(Means)]

  # Normalize weights
  sumWeights <- sum(Weights0)
  if (sumWeights != 1) {
    Weights0 <- Weights0 / sumWeights
  }

  # Calculate combined mean and SD
  Ns <- Weights0 * n
  MeanAll <- sum(Means0 * Weights0)

  # Calculate combined SD
  if (n == 1) {
    sc <- 0
  } else {
    qc <- sum((Ns - 1) * SDs0^2 + Ns * Means0^2)
    variance <- (qc - n * MeanAll^2) / (n - 1)
    sc <- sqrt(max(0, variance))  # Prevent negative variance from rounding errors
  }

  return(list(Mean = MeanAll, SD = sc))
}
