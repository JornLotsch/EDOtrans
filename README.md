# EDOtrans

Euclidean Distance-Optimized Data Transformation for Cluster Analysis

## Overview

The `EDOtrans` package implements the Euclidean distance-optimized (EDO) data transformation method for preprocessing data prior to cluster analysis. This transformation addresses limitations of standard z-standardization when working with complex, multimodal distributed variables.

## Background

The Euclidean metric is not scale invariant and can be inappropriate for complex distributed variables, negatively affecting cluster analysis results. The EDO transformation scales variables based on the standard deviations of dominant modes identified through Gaussian mixture modeling, optimizing distances for the Euclidean metric used in clustering algorithms.

Studies on artificial and biomedical datasets have shown that EDO transformation outperforms classical alternatives (untransformed data, z-transformation) in terms of cluster accuracy, adjusted Rand index, and Dunn's index.

## Installation

```r
# From CRAN
install.packages("EDOtrans")

# From GitHub
# install.packages("devtools")
devtools::install_github("JornLotsch/EDOtrans")
```

## Usage

```r
library(EDOtrans)

# Basic usage with automatic mode detection
result <- EDOtrans(Data = your_data)

# With predefined classes
result <- EDOtrans(Data = your_data, Cls = your_classes)

# Access transformed data
transformed_data <- result$DataEDO
scaling_factor <- result$EDOfactor
classes <- result$Cls
```

## Reference

Ultsch, A., & Lötsch, J. (2022). Euclidean distance-optimized data transformation for cluster analysis in biomedical data (EDOtrans). *BMC Bioinformatics*, 23(1), 233. https://doi.org/10.1186/s12859-022-04769-w

## License

License: GPL-3

## Related links

- [CRAN package page](https://cran.r-project.org/package=EDOtrans) 

- [imports cABC analysis R package](https://github.com/AndreHDev/cABC_Analysis) 
