#Calculates parameters of combined Gaussian modes

combinedModesParameters <- function(Means, SDs, Weights, n) {
  Means0 <- Means[order(Means)]
  SDs0 <- SDs[order(Means)]
  Weights0 <- Weights[order(Means)]
  sumWeights <- sum(Weights0)
  if (sumWeights != 1) {
    Weights0 <- Weights0 / sumWeights
  }
  Ns <- Weights0 * n
  MeanAll <- sum(Means0 * Weights0)
  qc <- sum((Ns - 1) * SDs0 ^ 2 + Ns * Means0 ^ 2)
  sc <- sqrt((qc - (n) * MeanAll ^ 2) / (n - 1))
  return(list(Mean = MeanAll, SD = sc))
}
