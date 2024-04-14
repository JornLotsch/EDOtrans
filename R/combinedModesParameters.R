# Calculates parameters of combined Gaussian modes
combinedModesParameters <- function(Means, SDs, Weights, n) {
  # Order the parameters by the means
  Means0 <- Means[order(Means)]
  SDs0 <- SDs[order(Means)]
  Weights0 <- Weights[order(Means)]

  # Normalize the weights if necessary
  sumWeights <- sum(Weights0)
  if (sumWeights != 1) {
    Weights0 <- Weights0 / sumWeights
  }

  # Calculate the number of samples per mode
  Ns <- Weights0 * n

  # Calculate the combined mean and standard deviation
  MeanAll <- sum(Means0 * Weights0)
  qc <- sum((Ns - 1) * SDs0^2 + Ns * Means0^2)
  sc <- sqrt((qc - (n) * MeanAll^2) / (n - 1))

  return(list(Mean = MeanAll, SD = sc))
}
