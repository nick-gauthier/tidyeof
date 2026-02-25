---
title: 'tidyeof: Empirical Orthogonal Functions and Statistical Downscaling for Spatiotemporal Data in R'
tags:
  - R
  - climate science
  - empirical orthogonal functions
  - canonical correlation analysis
  - downscaling
  - spatiotemporal data
authors:
  - name: Nicolas Gauthier
    orcid: 0000-0002-2225-5827
    affiliation: 1
affiliations:
  - name: Department of Anthropology & Florida Museum of Natural History, University of Florida, USA
    index: 1
date: 6 February 2026
bibliography: paper.bib
---

# Summary

`tidyeof` is an R package for Empirical Orthogonal Function (EOF) analysis and EOF-based statistical downscaling of spatiotemporal data. It provides a complete pipeline from raw climate fields to downscaled predictions: extracting dominant modes of variability via EOF decomposition, coupling predictor and response fields through Canonical Correlation Analysis (CCA), and validating results with temporal cross-validation. The package operates on `stars` [@pebesma2023] spatiotemporal arrays and integrates with the tidyverse ecosystem, making modern climate analysis methods accessible to researchers across disciplines.

# Statement of Need

EOF analysis is a foundational tool in climate science and related fields for decomposing spatiotemporal variability into interpretable modes [@hannachi2007]. When combined with CCA, it forms the basis of statistical downscaling — transferring information from coarse-resolution climate models or reanalyses to finer spatial scales [@vonStorch1999; @barnett1987]. These methods remain central to climate impact assessment, paleoclimate reconstruction, and regional climate analysis.

Despite the widespread use of these methods, no existing R package provides an integrated EOF-to-downscaling workflow built on modern spatial data infrastructure. In Python, `eofs` [@dawson2016] provides a mature EOF solver with NumPy, Iris, and xarray interfaces, and `xeofs` [@rieger2024] extends this with rotated EOF, Maximum Covariance Analysis, and dask support. However, neither includes a downscaling pipeline, and neither serves R users.

In R, existing options are fragmented. The `esd` package [@benestad2015] offers comprehensive statistical downscaling but predates `stars` and `sf` spatial objects and is not available on CRAN. The `wql` package provides a basic `eof()` function oriented toward water quality analysis. Packages in the `climate4R` bundle [@iturbide2019] include EOF-based downscaling but as part of a large, monolithic framework. Researchers commonly implement EOF-CCA pipelines from scratch using base R matrix operations, resulting in duplicated effort and inconsistent handling of spatial metadata, area weighting, and unit propagation.

`tidyeof` addresses this gap with a cohesive package that handles the full workflow — from climatology computation and anomalization through EOF extraction, CCA coupling, prediction, and cross-validation — while preserving spatial coordinates, units, and metadata throughout. By building on `stars` for spatiotemporal arrays and `sf` for geometry, the package integrates naturally with R's modern geospatial stack and the tidyverse grammar of data manipulation.

# Package Design

The package is organized around three S3 classes that mirror the analytical workflow:

- **`patterns`**: EOF spatial patterns, amplitude time series, eigenvalues, and climatology, computed via `patterns()` with options for area weighting, standardization, varimax rotation, and monthly anomalization. Large datasets trigger automatic truncated SVD via `irlba` [@baglama2005] for computational efficiency.

- **`coupled_patterns`**: CCA-coupled predictor and response EOF patterns, computed via `couple()`. The `predict()` method performs spatial downscaling by projecting new coarse-resolution data through the learned CCA transformation to produce fine-resolution fields.

- **`common_patterns`**: Joint EOF decomposition across multiple datasets via `common_patterns()`, ensuring consistent spatial patterns for cross-source comparison — particularly useful when coupling reanalyses with paleoclimate reconstructions or transferring models across climate products.

Each class supports standard R generics (`print`, `plot`, `summary`, `predict`, `[`, `$`) so that EOF objects behave like familiar R objects. Pattern visualization uses `ggplot2` with `scico` perceptually uniform color scales.

Diagnostics include scree plots with North et al. [-@north1982] error bars for mode separability, modified Rule N significance testing based on the Tracy-Widom distribution, teleconnection correlation maps with FDR-corrected significance contours, and canonical correlation visualization for coupled patterns.

Cross-validation functions (`tune_eof`, `tune_cca`) support data-driven selection of EOF truncation and CCA dimensionality using contiguous temporal folds, with pre-computed patterns enabling efficient grid search over hyperparameter combinations. This avoids the common pitfall of selecting modes by explained variance alone without assessing predictive skill [@wilks2011].

# Research Applications

An earlier version of the package was used in @gauthier2022 to downscale snowpack variability across the western United States, coupling large-scale atmospheric circulation patterns with station-level snow water equivalent observations via CCA. The current version supports several active research programs spanning climate science and archaeology. It provides the downscaling infrastructure for paleoclimate reconstruction in studies of medieval plague dynamics, coupling the PHYDA paleo-reanalysis [@steiger2018] with high-resolution CHELSA climatologies to drive epidemiological and demographic models at scales relevant to disease transmission. It is also used for temperature and precipitation downscaling over the US Southwest for archaeological demographic simulation. In each case, the package's cross-validation framework enables rigorous assessment of downscaling skill before using the predictions in downstream models.

# Acknowledgements

Development of `tidyeof` was supported by the Burroughs Wellcome Fund, NSF RISE award 2409068, and NSF award AGS-1803995.

# AI Usage Disclosure

Generative AI tools (Claude, Anthropic) were used during development for code review, test writing, and drafting portions of documentation. All AI-generated content was reviewed and edited by the author. The core algorithmic implementations were written by the author.

# References
