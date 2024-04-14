# Performs EDO transformation on given data and classes
#' @importFrom ABCanalysis ABCanalysis
#' @importFrom methods hasArg
#' @importFrom stats median sd runif
#' @importFrom opGMMassessment opGMMassessment
#' @export
EDOtrans <- function(Data, Cls, PlotIt = FALSE, FitAlg = "normalmixEM", Criterion = "LR",
                     MaxModes = 8, MaxCores = getOption("mc.cores", 2L), Seed) {
  # Check if data is provided
  if (!hasArg("Data")) {
    stop("EDOtrans: No data provided. Stopping.")
  }

  # Define a helper function to check if a vector is an empty integer vector
  is.integer0 <- function(x) {
    is.integer(x) && length(x) == 0L
  }

  # Set the seed if provided, otherwise use the current seed
  if (missing(Seed)) {
    ActualSeed <- as.integer(get_seed()[1])
  } else {
    ActualSeed <- Seed
  }

  # Determine the means, standard deviations, and weights based on the provided classes
  if (hasArg("Cls")) {
    if (length(Cls) != length(Data)) {
      stop("EDOtrans: Classes provided but unequal lengths of Data and Cls.")
    } else {
      Means0 <- tapply(Data, Cls, mean)
      SDs0 <- tapply(Data, Cls, sd)
      Weights0 <- tapply(Data, Cls, function(x) length(x) / length(Data))
    }
  } else {
    # Obtain classes via opGMMassessment
    warning("EDOtrans: Classes created using Gaussian mixture modeling.", call. = FALSE)
    GMMresults <- opGMMassessment(Data = Data, FitAlg = FitAlg, Criterion = Criterion,
                                  MaxModes = MaxModes, MaxCores = MaxCores, PlotIt = PlotIt, KS = FALSE,
                                  Seed = ActualSeed)
    Cls <- GMMresults$Cls
    Means0 <- GMMresults$Means
    SDs0 <- GMMresults$SDs
    Weights0 <- GMMresults$Weights
  }

  # Select the dominant groups
  if (length(Weights0) > 1) {
    WeightsABC <- ABCanalysis::ABCanalysis(as.vector(Weights0))
    if (is.integer0(WeightsABC$Aind) == FALSE) {
      Means0 <- Means0[WeightsABC$Aind]
      SDs0 <- SDs0[WeightsABC$Aind]
      Weights0 <- Weights0[WeightsABC$Aind]
      # Combine standard deviations from dominant groups
      nDom <- sum(Weights0) * length(Data)
      CombinedDominatGroupsParameters <- combinedModesParameters(Means = Means0,
                                                                 SDs = SDs0, Weights = Weights0, n = nDom)
      SDdomSq <- CombinedDominatGroupsParameters$SD * sqrt(2)
    } else {
      SDdomSq <- stats::median(SDs0) * sqrt(2)
    }
  } else {
    SDdomSq <- SDs0 * sqrt(2)
  }

  # Perform the EDO transformation
  DataEDOtrans <- Data / SDdomSq

  return(list(DataEDO = DataEDOtrans, EDOfactor = SDdomSq, Cls = Cls))
}
