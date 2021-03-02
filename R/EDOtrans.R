#Performs EDO transformation on given data and classes
#' @importFrom ABCanalysis ABCanalysis
#' @importFrom AdaptGauss BayesDecisionBoundaries
#' @importFrom methods hasArg
#' @importFrom stats dnorm median na.omit sd
#' @importFrom DistributionOptimization DistributionOptimization
#' @export
EDOtrans <- function(Data, Cls, Means, SDs, Weights, DO = FALSE, PlotGMM = FALSE) {
  if (hasArg("Data") == FALSE)
    stop("EDOtrans: No data provided. Stopping.")
  is.integer0 <- function(x)
  {
    is.integer(x) && length(x) == 0L
  }

  #Obtain standard deviations
  if (hasArg("Means") == TRUE & hasArg("SDs") == TRUE & hasArg("Weights") == TRUE) {
    if (length(c(Means, SDs, Weights)) %% 3 != 0)
      stop("EDOtrans: Unequal number of modes in parameters.")
    Means0 <- Means
    SDs0 <- SDs
    Weights0 <- Weights
    Boundaries <- c()
    Cls <- rep(1, length(Data))
    if (length(Means0) > 1) {
      Boundaries <- AdaptGauss::BayesDecisionBoundaries(Means = Means0, SDs = SDs0, Weights = Weights0)
      if (is.integer0(Boundaries) == FALSE)
        Boundaries <- Boundaries[Boundaries >= min(Means0) & Boundaries <= max(Means0)]
      if (length(Boundaries) > 0) {
        Cls <- cutGMM(x = Data, breaks = Boundaries)
      }
    }

  } else {
    if (hasArg("Cls") == TRUE) {
      if (length(Cls) != length(Data)) {
        stop("EDOtrans: Unequal legths of Data and Cls.")
      } else {
        Means0 <- tapply(X = Data, INDEX =  Cls, FUN = mean)
        SDs0 <- tapply(X = Data, INDEX =  Cls, FUN = sd)
        Weights0 <- tapply(X = Data, INDEX =  Cls, function(x) length(x) / length(Data))
      }
    } else {
      #Obtain classes via GMM
      warning("EDOtrans: Classes created using Gaussian mixtures.", call. = FALSE)
      if (hasArg("DO") == TRUE) DO = DO else DO = FALSE
      if (hasArg("PlotGMM") == TRUE) PlotGMM = PlotGMM else PlotGMM = FALSE
      GMMresults  <- GMMasessment(Data = Data, DO = DO, PlotIt = PlotGMM)
      Cls <- GMMresults$Cls
      Means0 <- GMMresults$Means
      SDs0 <- GMMresults$SDs
      Weights0 <- GMMresults$Weights
    }
  }

  #Selection of dominant groups
  if (length(Weights0) > 1) {
    WeightsABC <- ABCanalysis::ABCanalysis(as.vector(Weights0))
    if (is.integer0(WeightsABC$Aind) == FALSE) {
      Means0 <- Means0[WeightsABC$Aind]
      SDs0 <- SDs0[WeightsABC$Aind]
      Weights0 <- Weights0[WeightsABC$Aind]
      #Combine standard deviations  from dominant groups
      nDom <- sum(Weights0) * length(Data)
      CombinedDominatGroupsParameters <- combinedModesParameters(Means = Means0, SDs = SDs0, Weights = Weights0, n = nDom)
      SDdomSq <- CombinedDominatGroupsParameters$SD * sqrt(2)
    } else {
      SDdomSq <- median(SDs0) * sqrt(2)
    }
  } else {
    SDdomSq <- SDs0 * sqrt(2)
  }

  #Perform EDO transformation
  DataEDOtrans <- Data / SDdomSq

  return(list(DataEDO = DataEDOtrans, EDOfactor = SDdomSq, Cls = Cls))
}
