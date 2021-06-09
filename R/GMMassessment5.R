#Analysis of a Gaussian mixture structure in the data
#Statistical justification using likelihood ratio tests of GMM_M versus GMM_M-1, or AIC or BIC
#' @importFrom ClusterR GMM
#' @importFrom parallel detectCores
#' @importFrom pbmcapply pbmclapply
#' @importFrom AdaptGauss InformationCriteria4GMM LikelihoodRatio4Mixtures
#' @importFrom grDevices nclass.FD
#' @importFrom methods hasArg
#' @importFrom stats dnorm median na.omit sd ks.test
#' @importFrom DistributionOptimization DistributionOptimization
GMMasessment <-
  function(Data, DO = FALSE, PlotIt = FALSE, KS = FALSE, Criterion = "LR", MaxModes = 10, Seed) {
    if (!hasArg("Data"))
      stop("GMMasessment: No data.")
    if (length(Data) < 2)
      stop("GMMasessment: Too few data.")

    if (!missing(Seed)) {
      ActualSeed <- Seed
    } else {
      ActualSeed <- tail(get('.Random.seed', envir = globalenv()), 1)
    }

    is.integer0 <- function(x) {
      is.integer(x) && length(x) == 0L
    }

    idBestGMM_AICBIC <- function(GMMdata, GMMfit, Criterion) {
      AICBIC <- lapply(1:length(GMMfit), function(x) {
        AICBICi <- AdaptGauss::InformationCriteria4GMM(
          Data = GMMdata,
          Means = lapply(GMMfit, "[[", 2)[[x]][, 1],
          SDs = lapply(GMMfit, "[[", 2)[[x]][, 2],
          Weights = lapply(GMMfit, "[[", 2)[[x]][, 3]
        )
        return(AICBICi)
      })
      AICBIC <- unlist(lapply(1:length(GMMfit), function(x) {
        switch(Criterion,
               AIC = { AICBIC <- AICBIC[[x]]$AIC },
               BIC = { AICBIC <- AICBIC[[x]]$BIC })
        return(AICBIC)
      }))
      firstBestGMM <- 1
      for (i in 2:length(GMMfit)) {
        if (AICBIC[i] < AICBIC[i - 1])
          firstBestGMM <- i
        else
          break
      }
      minBestGMM <- which.min(AICBIC)
      if (firstBestGMM != minBestGMM) {
        require("twosamples")
        set.seed(ActualSeed)
        Pred <-
          CreateGMM(
            Means = lapply(GMMfit, "[[", 2)[[firstBestGMM]][, 1],
            SDs = lapply(GMMfit, "[[", 2)[[firstBestGMM]][, 2],
            Weights = lapply(GMMfit, "[[", 2)[[firstBestGMM]][, 3],
            n = 1000
          )$Data
        Pred <- Pred[Pred >= min(GMMdata) & Pred <= max(GMMdata)]
        KSfirst <- ks.test(x = GMMdata, y = Pred)$statistic

        set.seed(ActualSeed)
        Pred <-
          CreateGMM(
            Means = lapply(GMMfit, "[[", 2)[[minBestGMM]][, 1],
            SDs = lapply(GMMfit, "[[", 2)[[minBestGMM]][, 2],
            Weights = lapply(GMMfit, "[[", 2)[[minBestGMM]][, 3],
            n = 1000
          )$Data
        Pred <- Pred[Pred >= min(GMMdata) & Pred <= max(GMMdata)]
        KSmin <- ks.test(x = GMMdata, y = Pred)$statistic

        if (KSfirst < KSmin)
          BestGMM <- firstBestGMM
        else
          BestGMM <- minBestGMM
      } else
        BestGMM <- firstBestGMM
      return(BestGMM)
    }


    idBestGMM_LR <- function(GMMdata, GMMfit) {
      LRi <- c(1, unlist(lapply(2:MaxModes, function(x) {
        AdaptGauss::LikelihoodRatio4Mixtures(
          Data = GMMdata,
          NullMixture = lapply(GMMfit, "[[", 2)[[x - 1]],
          OneMixture = lapply(GMMfit, "[[", 2)[[x]],
          PlotIt = FALSE
        )$Pvalue
      })))

      LR1 <- c(1, unlist(lapply(2:MaxModes, function(x) {
        AdaptGauss::LikelihoodRatio4Mixtures(
          Data = GMMdata,
          NullMixture = lapply(GMMfit, "[[", 2)[[1]],
          OneMixture = lapply(GMMfit, "[[", 2)[[x]],
          PlotIt = FALSE
        )$Pvalue
      })))

      firstBestGMM <- 1
      for (i in 2:length(GMMfit)) {
        if (LRi[i] < 0.05)
          firstBestGMM <- i
        else
          break
      }
      minBestGMM <- which.min(LR1)
      if (LR1[minBestGMM] >= 0.05) minBestGMM <- firstBestGMM

      if (firstBestGMM != minBestGMM) {
        set.seed(ActualSeed)
        Pred <-
          CreateGMM(
            Means = lapply(GMMfit, "[[", 2)[[firstBestGMM]][, 1],
            SDs = lapply(GMMfit, "[[", 2)[[firstBestGMM]][, 2],
            Weights = lapply(GMMfit, "[[", 2)[[firstBestGMM]][, 3],
            n = 1000
          )$Data
        Pred <- Pred[Pred >= min(GMMdata) & Pred <= max(GMMdata)]
        KSfirst <- ks.test(x = GMMdata, y = Pred)$statistic

        set.seed(ActualSeed)
        Pred <-
          CreateGMM(
            Means = lapply(GMMfit, "[[", 2)[[minBestGMM]][, 1],
            SDs = lapply(GMMfit, "[[", 2)[[minBestGMM]][, 2],
            Weights = lapply(GMMfit, "[[", 2)[[minBestGMM]][, 3],
            n = 1000
          )$Data
        Pred <- Pred[Pred >= min(GMMdata) & Pred <= max(GMMdata)]
        KSmin <- ks.test(x = GMMdata, y = Pred)$statistic

        if (KSfirst < KSmin)
          BestGMM <- firstBestGMM
        else
          BestGMM <- minBestGMM
      } else
        BestGMM <- firstBestGMM
      return(BestGMM)
    }


    GMMdata <- Data
    MaxModes <- MaxModes
    list.of.Modes <- 1:MaxModes

    #Do the GMM fit
    if (.Platform$OS.type != "windows" &
        MaxModes > 1 & DO == TRUE) {
      require("parallel")
      require("pbmcapply")
      chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
      if (nzchar(chk) && chk == "TRUE") {
        num_workers <- 2L
      } else {
        num_workers <- detectCores()
      }
      nProc <- min(num_workers - 1, MaxModes)

      #GMM fit using a genetic algorithm
      GMMfit <- pbmclapply(list.of.Modes, function(x) {
        GMMfit_Mode <-
          DistributionOptimization::DistributionOptimization(
            Data = GMMdata,
            Modes = list.of.Modes[x],
            Monitor = 0,
            CrossoverRate = .9,
            ErrorMethod = "chisquare",
            Seed = ActualSeed
          )
        Mixtures <-
          cbind(GMMfit_Mode$Means,
                GMMfit_Mode$SDs,
                GMMfit_Mode$Weights)
        return(list(GMMfit_Mode, Mixtures))
      }, mc.cores = nProc)
    } else {
      if (DO == FALSE) {
        #GMM fit using EM
        GMMfit <- lapply(list.of.Modes, function(x) {
          GMMfit_Mode <-
            ClusterR::GMM(
              data = data.frame(GMMdata),
              gaussian_comps = list.of.Modes[x],
              dist_mode = "eucl_dist"
            )
          Mixtures <-
            cbind(
              GMMfit_Mode$centroids,
              sqrt(GMMfit_Mode$covariance_matrices),
              GMMfit_Mode$weights
            )
          return(list(GMMfit_Mode, Mixtures))
        })
      } else {
        #GMM fit using a genetic algorithm
        GMMfit <- lapply(list.of.Modes, function(x) {
          GMMfit_Mode <-
            DistributionOptimization::DistributionOptimization(
              Data = GMMdata,
              Modes = list.of.Modes[x],
              Monitor = 0,
              CrossoverRate = .9,
              ErrorMethod = "chisquare",
              Seed = ActualSeed
            )
          Mixtures <-
            cbind(GMMfit_Mode$Means,
                  GMMfit_Mode$SDs,
                  GMMfit_Mode$Weights)
          return(list(GMMfit_Mode, Mixtures))
        })
      }
    }

    #Identify best fit based on selected criterion
    switch(Criterion,
           BIC = { BestGMM <- idBestGMM_AICBIC(GMMdata, GMMfit, Criterion) },
           AIC = { BestGMM <- idBestGMM_AICBIC(GMMdata, GMMfit, Criterion) },
           LR = { BestGMM <- idBestGMM_LR(GMMdata, GMMfit) })

    if (DO == FALSE) {
      Means <- as.vector(GMMfit[[BestGMM]][[1]]$centroids)
      SDs <- sqrt(as.vector(GMMfit[[BestGMM]][[1]]$covariance_matrices))
      Weights <- as.vector(GMMfit[[BestGMM]][[1]]$weights)
    } else {
      Means <- as.vector(GMMfit[[BestGMM]][[1]]$Means)
      SDs <- as.vector(GMMfit[[BestGMM]][[1]]$SDs)
      Weights <- as.vector(GMMfit[[BestGMM]][[1]]$Weights)
    }

    #Calculate Bayes boundaries
    Boundaries <- c()
    Classes <- rep(1, length(GMMdata))
    if (BestGMM > 1) {
      Boundaries <-
        AdaptGauss::BayesDecisionBoundaries(Means = Means,
                                            SDs = SDs,
                                            Weights = Weights)
      if (is.integer0(Boundaries) == FALSE)
        Boundaries <-
          Boundaries[Boundaries >= min(Means) & Boundaries <= max(Means)]
      if (length(Boundaries) > 0)
        Classes <- cutGMM(x = GMMdata, breaks = Boundaries)
    }

    #Do Kolmogorov-Smirnov test
    if (KS == TRUE) {
      set.seed(ActualSeed)
      Pred <-
        CreateGMM(
          Means = lapply(GMMfit, "[[", 2)[[BestGMM]][, 1],
          SDs = lapply(GMMfit, "[[", 2)[[BestGMM]][, 2],
          Weights = lapply(GMMfit, "[[", 2)[[BestGMM]][, 3],
          n = 1000
        )$Data
      Pred <- Pred[Pred >= min(GMMdata) & Pred <= max(GMMdata)]
      KStest <- ks.test(x = GMMdata, y = Pred)
    } else
      KStest <- NA

    #Prepare plot
    p1 <-
      GMMplotGG(Data = GMMdata, Means = Means, SDs = SDs, Weights = Weights, Hist = TRUE)
    if (PlotIt == TRUE)
      print(p1)
    return(
      list(
        Cls = Classes,
        Means = Means,
        SDs = SDs,
        Weights = Weights,
        Boundaries = Boundaries,
        Plot = p1,
        KS = KStest
      )
    )
  }
