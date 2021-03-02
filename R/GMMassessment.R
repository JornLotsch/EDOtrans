#Analysis of a Gaussian mixture structure in the data
#Statistical justification using likelihood ratio tests of GMM_M versus GMM_M-1
#' @importFrom ClusterR GMM
#' @importFrom AdaptGauss InformationCriteria4GMM
#' @importFrom grDevices nclass.FD
#' @importFrom methods hasArg
#' @importFrom stats dnorm median na.omit sd
#' @importFrom DistributionOptimization DistributionOptimization
GMMasessment <- function(Data, DO = FALSE, PlotIt = FALSE) {
  if (!hasArg("Data"))
    stop("GMMasessment: No data.")
  if (length(Data) < 2)
    stop("GMMasessment: Too few data.")
  is.integer0 <- function(x)
  {
    is.integer(x) && length(x) == 0L
  }

  GMMdata <- Data
  MaxModes <- 10

  if (DO == FALSE) {
    #GMM fit using EM
    list.of.Modes <- 1:MaxModes
    GMMfit <- lapply(1:MaxModes, function(x) {
      GMMfit_Mode <- ClusterR::GMM(data = data.frame(GMMdata), gaussian_comps = list.of.Modes[x], dist_mode = "eucl_dist")
      GMMfit_Mode_FitGuete <- AdaptGauss::InformationCriteria4GMM(GMMdata,
                                                                  Means = GMMfit_Mode$centroids,
                                                                  SDs = sqrt(GMMfit_Mode$covariance_matrices),
                                                                  Weights = GMMfit_Mode$weights)$BIC
      return(list(GMMfit_Mode, GMMfit_Mode_FitGuete))
    })
    BIC <- unlist(lapply(GMMfit, "[[", 2))
    BestGMM <- 1
    for (i in 2:MaxModes) {
      if (BIC[i] < BIC[i - 1])
        BestGMM <- i
      else
        break
    }
    Means = as.vector(GMMfit[[BestGMM]][[1]]$centroids)
    SDs = sqrt(as.vector(GMMfit[[BestGMM]][[1]]$covariance_matrices))
    Weights = as.vector(GMMfit[[BestGMM]][[1]]$weights)
  } else {
    #GMM fit using DistributionOptimization
    BICGMM <- vector()
    for (i in 1:MaxModes)
    {
      assign(
        paste0("GMMfit_Modes", i),
        DistributionOptimization::DistributionOptimization(Data = GMMdata, Modes = i, Monitor = 0, ErrorMethod = "chisquare"))
      assign(
        paste0("GMMfit_Modes", i, "_FitGuete"),
        AdaptGauss::InformationCriteria4GMM(Data = GMMdata,
                                            Means = get(paste0("GMMfit_Modes", i))$Means,
                                            SDs = get(paste0("GMMfit_Modes", i))$SDs,
                                            Weights = get(paste0("GMMfit_Modes", i))$Weights))
      if (i == 1) {
        assign("BICGMM", get(paste0("GMMfit_Modes", i, "_FitGuete"))$BIC)
        BestGMM <- 1
      } else {
        if (get(paste0("GMMfit_Modes", i, "_FitGuete"))$BIC < BICGMM) {
          BICGMM <- get(paste0("GMMfit_Modes", i, "_FitGuete"))$BIC
          BestGMM <- i
        } else
          break
      }
    }
    Means = get(paste0("GMMfit_Modes", BestGMM))$Means
    SDs = get(paste0("GMMfit_Modes", BestGMM))$SDs
    Weights = get(paste0("GMMfit_Modes", BestGMM))$Weights
  }
  Boundaries <- c()
  Classes <- rep(1, length(GMMdata))
  if (BestGMM > 1) {
    Boundaries <- AdaptGauss::BayesDecisionBoundaries(Means = Means, SDs = SDs, Weights = Weights)
    if (is.integer0(Boundaries) == FALSE)
      Boundaries <- Boundaries[Boundaries >= min(Means) & Boundaries <= max(Means)]
    if (length(Boundaries) > 0) {
      Classes <- cutGMM(x = GMMdata, breaks = Boundaries)
    }
  }
  p1 <- GMMplotGG(Data = Data, Means = Means, SDs = SDs, Weights = Weights, Hist = TRUE)
  if(PlotIt == TRUE) print(p1)
  return(list(Cls = Classes, Means = Means, SDs = SDs, Weights = Weights, Boundaries = Boundaries, Plot = p1))
}
