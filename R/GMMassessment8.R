#Analysis of a Gaussian mixture structure in the data
#Statistical justification using likelihood ratio tests
#of GMM_M versus GMM_M-1, or AIC or BIC
#' @importFrom ClusterR GMM
#' @importFrom grDevices nclass.FD
#' @importFrom methods hasArg
#' @importFrom utils tail
#' @importFrom cluster clusGap
#' @importFrom stats ks.test
#' @importFrom DistributionOptimization DistributionOptimization
GMMassessment <-
  function(Data, DO = FALSE, PlotIt = FALSE, KS = FALSE, Criterion = "BIC",
           Razor = TRUE, MaxModes = 10, MaxCores = 28, Seed) {
    if (!hasArg("Data"))
      stop("GMMassessment: No data.")
    if (length(Data) < 2)
      stop("GMMassessment: Too few data.")

    DataOrigLength <- length(as.vector(Data))
    DataOrignoNA <- which(!is.na(Data) & !is.infinite(Data))
    Data <- as.vector(Data[DataOrignoNA])
    n <- ifelse(length(Data) < 1000, 1000, 2 * length(Data))

    if (!missing(Seed)) {
      ActualSeed <- Seed
    } else {
      ActualSeed <- tail(get(".Random.seed", envir = globalenv()), 1)
    }

    is.integer0 <- function(x) {
      is.integer(x) && length(x) == 0L
    }

    pam1 <- function(x,k) {
      list(cluster = pam(x, 3, metric = "euclidean", stand = FALSE, cluster.only=TRUE))
    }

    GMMdata <- Data
    MaxModes <- MaxModes

    gsPam1 <- clusGap(cbind(GMMdata,GMMdata), FUN = pam1, K.max = MaxModes, B = 60, verbose = F, spaceH0 = "original")
    M <- which.max(gsPam1$Tab[,4])

    #Do the GMM fit
    if (DO == FALSE) {
      #GMM fit using EM
      GMMfit <- ClusterR::GMM(
        data = data.frame(GMMdata),
        gaussian_comps = M,
        dist_mode = "eucl_dist"
      )
      Means <- GMMfit$centroids
      SDs <- sqrt(GMMfit$covariance_matrices)
      Weights <- GMMfit$weights
    } else {
      #GMM fit using a genetic algorithm
      GMMfit <- DistributionOptimization::DistributionOptimization(
        Data = GMMdata,
        Modes = M,
        Monitor = 0,
        CrossoverRate = .9,
        ErrorMethod = "chisquare",
        Seed = ActualSeed
      )
      Means <- GMMfit$Means
      SDs <-GMMfit$SDs
      Weights <- GMMfit$Weights
    }

    #Calculate Bayes boundaries
    Boundaries <- c()
    ClassesB <- rep(1, length(GMMdata))
    if (M > 1) {
      Boundaries <-
        AdaptGauss::BayesDecisionBoundaries(Means = Means,
                                            SDs = SDs,
                                            Weights = Weights)
      if (is.integer0(Boundaries) == FALSE)
        Boundaries <-
          Boundaries[Boundaries >= min(Means) & Boundaries <= max(Means)]
      if (length(Boundaries) > 0)
        ClassesB <- cutGMM(x = GMMdata, breaks = Boundaries)
    }

    Classes <- rep(NA, length(DataOrigLength))
    Classes[DataOrignoNA] <- ClassesB

    #Do Kolmogorov-Smirnov test
    if (KS == TRUE) {
      set.seed(ActualSeed)
      Pred <-
        CreateGMM(
          Means = Means, SDs = SDs, Weights = Weights,
          n = n
        )$Data
      KStest <- suppressWarnings(ks.test(x = GMMdata, y = Pred))
    } else
      KStest <- NA

    #Prepare plot
    p1 <-
      GMMplotGG(Data = GMMdata, Means = Means, SDs = SDs,
                Weights = Weights, Hist = TRUE)
    if (PlotIt == TRUE)
      print(p1)

    #Return results
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

