# Performs EDO transformation on given data and classes
#' @importFrom ABCanalysis ABCanalysis
#' @importFrom methods hasArg
#' @importFrom stats median sd runif
#' @importFrom opGMMassessment opGMMassessment
#' @export
EDOtrans <- function(Data, Cls, PlotIt = FALSE, FitAlg = "normalmixEM", Criterion = "LR",
                     MaxModes = 8, MaxCores = getOption("mc.cores", 2L), Seed = "simple") {

  # Check Data input
  if (missing(Data)) stop("EDOtrans: No data provided. Stopping.")
  if (length(Data) == 0) stop("EDOtrans: Data is empty.")
  if (is.factor(Data)) Data <- as.character(Data)
  if (is.numeric(Data)) {
    Data <- as.numeric(Data)
  } else if (is.integer(Data)) {
    Data <- as.numeric(Data)
  } else if (is.character(Data) || is.factor(Data)) {
    x <- as.character(Data)
    ok <- !is.na(suppressWarnings(as.numeric(x))) | is.na(x)
    if (!all(ok)) stop("EDOtrans: Data must be numeric or coercible to numeric.")
    Data <- as.numeric(x)
  } else {
    stop("EDOtrans: Data must be numeric or coercible to numeric.")
  }
  if (all(is.na(Data))) stop("EDOtrans: Data contains only NA values.")

  # Check other input
  if (MaxModes < 1) {
    warning("EDOtrans: MaxModes < 1, setting to 1.", call. = FALSE)
    MaxModes <- 1
  }
  if (!FitAlg %in% c("ClusterRGMM", "densityMclust", "DO", "MCMC", "normalmixEM")) {
    warning("EDOtrans: Invalid FitAlg '", FitAlg, "', using 'normalmixEM'.", call. = FALSE)
    FitAlg <- "normalmixEM"
  }
  if (!Criterion %in%  c("AIC", "BIC", "FM", "GAP", "LR", "NbClust", "SI")) {
    warning("EDOtrans: Invalid Criterion '", Criterion, "', using 'LR'.", call. = FALSE)
    Criterion <- "LR"
  }

  # Seed handling with three clear options
  # Set the seed if provided, otherwise use the current seed
  if (missing(Seed)) {
    ActualSeed <- as.integer(get_seed())
  } else {
    if (is.numeric(Seed)) {
      # Option 1: Use provided integer seed directly
      ActualSeed <- as.integer(Seed)
      set.seed(ActualSeed)
    } else if (is.character(Seed)) {
      ActualSeed <- switch(Seed,
                           "auto" = as.integer(get_seed()),           # Complex seed recovery
                           "simple" = {                               # Simple reproducible seed (default)
                             temp_seed <- sample(1:100, 1)
                             warning(paste0("EDOtrans: Seed set at ", temp_seed, "."), call. = FALSE)
                             temp_seed
                           },
                           stop("Invalid Seed input. Use 'auto', 'simple', or an integer.")
      )
    } else {
      # Fallback for backward compatibility
      ActualSeed <- as.integer(get_seed())
    }
  }

  # Main part If classes are specified, transformation is done based on the
  # classes, otherwise the modality is checked automatically.
  if  (!missing(Cls) && !is.null(Cls)) {
    if (length(Cls) != length(Data)) {
      stop("EDOtrans: Classes provided but unequal lengths of Data and Cls.")
    } else {
      Means0 <- tapply(X = Data, INDEX = Cls, function(x) mean(x, na.rm = TRUE))
      SDs0 <- tapply(X = Data, INDEX = Cls, function(x) sd(x, na.rm = TRUE))
      Weights0 <- tapply(X = Data, INDEX = Cls, function(x) length(x)/length(Data))
      # Check for empty classes
      if (any(is.na(Means0)) || any(is.na(SDs0))) {
        stop("EDOtrans: Some classes have no data or only NA values.")
      }
      if (any(Weights0 < 0) || any(is.na(Weights0))) {
        stop("EDOtrans: Invalid weights computed from classes.")
      }
    }
  } else {
    if (MaxModes == 1) {
      Cls <- rep(1, length(Data))
      Means0 <- mean(Data, na.rm = TRUE)
      SDs0 <- sd(Data, na.rm = TRUE)
      Weights0 <- 1
    } else {
      # Obtain classes via opGMMassessment
      warning("EDOtrans: Classes created using Gaussian mixture modeling.", call. = FALSE)
      GMMresults <- opGMMassessment::opGMMassessment(Data = Data, FitAlg = FitAlg, Criterion = Criterion,
                                                     MaxModes = MaxModes, MaxCores = MaxCores, PlotIt = PlotIt, KS = FALSE,
                                                     Seed = ActualSeed)
      Cls <- GMMresults$Cls
      Means0 <- GMMresults$Means
      SDs0 <- GMMresults$SDs
      Weights0 <- GMMresults$Weights
    }
  }

  # Selection of dominant groups
  is.integer0 <- function(x) {
    is.integer(x) && length(x) == 0L
  }

  # Validate Weights0 before ABC analysis
  if (any(is.na(Weights0)) || any(Weights0 < 0)) {
    stop("EDOtrans: Invalid weights (NA or negative values).")
  }

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

  # Perform EDO transformation
  if (is.na(SDdomSq) || SDdomSq <= 0) {
    stop("EDOtrans: invalid scaling factor computed.")
  } else {
    DataEDOtrans <- Data/SDdomSq
  }

  return(list(DataEDO = DataEDOtrans, EDOfactor = SDdomSq, Cls = Cls))
}
