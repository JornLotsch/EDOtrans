# Performs EDO transformation on given data and classes
#' @importFrom ABCanalysis ABCanalysis
#' @importFrom methods hasArg
#' @importFrom stats median sd
#' @importFrom opGMMassessment opGMMassessment
#' @export
EDOtrans <- function(Data, Cls, PlotIt = FALSE, FitAlg = "normalmixEM", Criterion = "LR",
  MaxModes = 8, MaxCores = getOption("mc.cores", 2L), Seed) {
  if (hasArg("Data") == FALSE) {
    stop("EDOtrans: No data provided. Stopping.")
  }

  is.integer0 <- function(x) {
    is.integer(x) && length(x) == 0L
  }

  if (!missing(Seed)) {
    ActualSeed <- Seed
  } else {
    ActualSeed <- tail(get(".Random.seed", envir = globalenv()), 1)
  }

  # Main part If classes are specified, transformation is done based on the
  # classes, otherwise the modality is checked automatically.
  if (hasArg("Cls") == TRUE) {
    if (length(Cls) != length(Data)) {
      stop("EDOtrans: Classes provided but unequal legths of Data and Cls.")
    } else {
      Means0 <- tapply(X = Data, INDEX = Cls, FUN = mean)
      SDs0 <- tapply(X = Data, INDEX = Cls, FUN = sd)
      Weights0 <- tapply(X = Data, INDEX = Cls, function(x) length(x)/length(Data))
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

  # Selection of dominant groups
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
      SDdomSq <- median(SDs0) * sqrt(2)
    }
  } else {
    SDdomSq <- SDs0 * sqrt(2)
  }

  # Perform EDO transformation
  DataEDOtrans <- Data/SDdomSq

  return(list(DataEDO = DataEDOtrans, EDOfactor = SDdomSq, Cls = Cls))
}
