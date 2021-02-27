#Calculates parameters of combined Gaussian modes
#Modifed to > 2 modes from
#https://math.stackexchange.com/questions/2971315/how-do-i-combine-standard-deviations-of-two-groups
#@MISC {2971522,
#TITLE = {How do I combine standard deviations of two groups?},
#AUTHOR = {BruceET (https://math.stackexchange.com/users/221800/bruceet)},
#HOWPUBLISHED = {Mathematics Stack Exchange},
#NOTE = {URL:https://math.stackexchange.com/q/2971522 (version: 2018-10-26)},
#EPRINT = {https://math.stackexchange.com/q/2971522},
#URL = {https://math.stackexchange.com/q/2971522}
#}

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
