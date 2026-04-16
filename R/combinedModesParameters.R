# Calculates parameters of combined Gaussian modes
combinedModesParameters <- function(Means, SDs, Weights, n) {
<<<<<<< HEAD
  # Order the parameters by the means
=======
  # Validation
  if (length(Means) != length(SDs) || length(Means) != length(Weights)) {
    stop("Means, SDs, and Weights must have the same length")
  }
  if (length(n) != 1 || !is.numeric(n) || n <= 0) {
    stop("n must be a positive number")
  }

  # Order by means
>>>>>>> 8c002aa (get_seed new, error handling improved)
  Means0 <- Means[order(Means)]
  SDs0 <- SDs[order(Means)]
  Weights0 <- Weights[order(Means)]

<<<<<<< HEAD
  # Normalize the weights if necessary
=======
  # Normalize weights
>>>>>>> 8c002aa (get_seed new, error handling improved)
  sumWeights <- sum(Weights0)
  if (sumWeights != 1) {
    Weights0 <- Weights0 / sumWeights
  }

<<<<<<< HEAD
  # Calculate the number of samples per mode
  Ns <- Weights0 * n

  # Calculate the combined mean and standard deviation
  MeanAll <- sum(Means0 * Weights0)
  qc <- sum((Ns - 1) * SDs0^2 + Ns * Means0^2)
  sc <- sqrt((qc - (n) * MeanAll^2) / (n - 1))
=======
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
>>>>>>> 8c002aa (get_seed new, error handling improved)

  return(list(Mean = MeanAll, SD = sc))
}
