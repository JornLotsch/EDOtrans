#Analysis of a Gaussian mixture structure in the data
#Statistical justification using likelihood ratio tests of GMM_M versus GMM_M-1
#' @importFrom ClusterR GMM
#' @importFrom AdaptGauss InformationCriteria4GMM LikelihoodRatio4Mixtures
#' @importFrom grDevices nclass.FD
#' @importFrom methods hasArg
#' @importFrom stats dnorm median na.omit sd
#' @importFrom DistributionOptimization DistributionOptimization
GMMasessment <- function(Data, DO = FALSE, PlotIt = FALSE, Criterion = "BIC", MaxModes = 10) {
  if (!hasArg("Data"))
    stop("GMMasessment: No data.")
  if (length(Data) < 2)
    stop("GMMasessment: Too few data.")
  is.integer0 <- function(x) {
    is.integer(x) && length(x) == 0L
  }

  GMMdata <- Data
  MaxModes <- MaxModes
  list.of.Modes <- 1:MaxModes

  #Do the GMM fit
  if (DO == FALSE) {
    #GMM fit using EM
    GMMfit <- lapply(list.of.Modes, function(x) {
      GMMfit_Mode <- ClusterR::GMM(data = data.frame(GMMdata), gaussian_comps = list.of.Modes[x], dist_mode = "eucl_dist")
      Mixtures = cbind(GMMfit_Mode$centroids, sqrt(GMMfit_Mode$covariance_matrices), GMMfit_Mode$weights)
      return(list(GMMfit_Mode, Mixtures))
    })

    #Identify best fit based on selected criterion
    switch(Criterion,
            BIC = {
      BIC <- unlist(lapply(list.of.Modes, function(x) {
        BICi <- AdaptGauss::InformationCriteria4GMM(
                  Data = GMMdata,
                  Means = lapply(GMMfit, "[[", 2)[[x]][, 1],
                  SDs = lapply(GMMfit, "[[", 2)[[x]][, 2],
                  Weights = lapply(GMMfit, "[[", 2)[[x]][, 3]
                )$BIC
        return(BICi)
      }))
      BestGMM <- 1
      for (i in 2:MaxModes) {
        if (BIC[i] < BIC[i - 1])
          BestGMM <- i
        else
          break
      }
    },
            AIC = {
      AIC <- unlist(lapply(list.of.Modes, function(x) {
        AICi <- AdaptGauss::InformationCriteria4GMM(
                  Data = GMMdata,
                  Means = lapply(GMMfit, "[[", 2)[[x]][, 1],
                  SDs = lapply(GMMfit, "[[", 2)[[x]][, 2],
                  Weights = lapply(GMMfit, "[[", 2)[[x]][, 3]
                )$AIC
        return(AICi)
      }))
      BestGMM <- 1
      for (i in 2:MaxModes) {
        if (AIC[i] < AIC[i - 1])
          BestGMM <- i
        else
          break
      }
    },
            LR = {
      BestGMM <- 1
      for (i in 1:(MaxModes - 1)) {
        LRp <- AdaptGauss::LikelihoodRatio4Mixtures(
                  Data = GMMdata,
                  NullMixture = lapply(GMMfit, "[[", 2)[[i]],
                  OneMixture = lapply(GMMfit, "[[", 2)[[i + 1]],
                  PlotIt = FALSE
                )$Pvalue
        if (LRp < 0.05)
          BestGMM <- i
        else
          break
      }
    })

    Means = as.vector(GMMfit[[BestGMM]][[1]]$centroids)
    SDs = sqrt(as.vector(GMMfit[[BestGMM]][[1]]$covariance_matrices))
    Weights = as.vector(GMMfit[[BestGMM]][[1]]$weights)
  } else {
    #GMM fit using DistributionOptimization, for speed, as a loop that breaks when no further improvement
    for (i in 1:MaxModes) {
      assign(
        paste0("GMMfit_Modes", i),
        DistributionOptimization::DistributionOptimization(
          Data = GMMdata,
          Modes = i,
          Monitor = 0,
          ErrorMethod = "chisquare"
        )
      )

      #Identify best fit based on selected criterion
      assign(
        paste0("GMMfit_Modes", i, "_FitGuete"),
        AdaptGauss::InformationCriteria4GMM(
          Data = GMMdata,
          Means = get(paste0("GMMfit_Modes", i))$Means,
          SDs = get(paste0("GMMfit_Modes", i))$SDs,
          Weights = get(paste0("GMMfit_Modes", i))$Weights
        )
      )
      if (i == 1) {
        assign("BICGMM", get(paste0("GMMfit_Modes", i, "_FitGuete"))$BIC)
        assign("AICGMM", get(paste0("GMMfit_Modes", i, "_FitGuete"))$AIC)
        BestGMM <- 1
      } else {
        switch(Criterion,
                BIC = {
          if (get(paste0("GMMfit_Modes", i, "_FitGuete"))$BIC < BICGMM) {
            BICGMM <- get(paste0("GMMfit_Modes", i, "_FitGuete"))$BIC
            BestGMM <- i
          } else
            break
        },
                AIC = {
        if (get(paste0("GMMfit_Modes", i, "_FitGuete"))$AIC < AICGMM) {
          AICGMM <- get(paste0("GMMfit_Modes", i, "_FitGuete"))$AIC
          BestGMM <- i
        } else
          break
      },
                LR = {
      NullMixture = cbind(get(paste0("GMMfit_Modes", i - 1))$Means,
                                      get(paste0("GMMfit_Modes", i - 1))$SDs,
                                      get(paste0("GMMfit_Modes", i - 1))$Weights)
      OneMixture = cbind(get(paste0("GMMfit_Modes", i))$Means,
                                     get(paste0("GMMfit_Modes", i))$SDs,
                                     get(paste0("GMMfit_Modes", i))$Weights)
      LRp <-
                    AdaptGauss::LikelihoodRatio4Mixtures(
                      Data = GMMdata,
                      NullMixture = NullMixture,
                      OneMixture = OneMixture,
                      PlotIt = FALSE
                    )
      if (LRp < 0.05)
        BestGMM <- i
      else
        break
    })
      }
      Means = get(paste0("GMMfit_Modes", BestGMM))$Means
      SDs = get(paste0("GMMfit_Modes", BestGMM))$SDs
      Weights = get(paste0("GMMfit_Modes", BestGMM))$Weights
    }
  }

  #Calculate Bayes boundaries
  Boundaries <- c()
  Classes <- rep(1, length(GMMdata))
  if (BestGMM > 1) {
    Boundaries <- AdaptGauss::BayesDecisionBoundaries(Means = Means,
                                                        SDs = SDs,
                                                        Weights = Weights)
    if (is.integer0(Boundaries) == FALSE)
      Boundaries <- Boundaries[Boundaries >= min(Means) & Boundaries <= max(Means)]
    if (length(Boundaries) > 0) {
      Classes <- cutGMM(x = GMMdata, breaks = Boundaries)
    }
  }

  #Prepare plot
  p1 <- GMMplotGG(
      Data = Data,
      Means = Means,
      SDs = SDs,
      Weights = Weights,
      Hist = TRUE
    )
  if (PlotIt == TRUE)
    print(p1)
  return(
      list(
        Cls = Classes,
        Means = Means,
        SDs = SDs,
        Weights = Weights,
        Boundaries = Boundaries,
        Plot = p1
      )
    )
}
