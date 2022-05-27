#!/usr/bin/env Rscript

# Ensure proper setup
installed <- installed.packages()[, 1]
# as of 2022-05-27, this is needed for patched version of
# readBrukerFlexData::readBrukerFlexData

if (! "readBrukerFlexData" %in% installed) {
  cat("Seting up readBruckerFlexData\n")
  if (!"remotes" %in% installed) {
    install.packages("remotes", repos = "https://cran.r-project.org")
  }
  remotes::install_github("sgibb/readBrukerFlexData@issue3-ctof2calibration")
}
library("readBrukerFlexData")

for (pkg in c("MALDIquant",
              "MALDIquantForeign",
              "optparse",
              "magrittr",
              "readr",
              "purrr")) {
  if (!pkg %in% installed) {
    cat("Installing", pkg, "\n")
    install.packages(pkg, repos = "https://cran.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}


arguments <- function() {
  parser <-
    OptionParser() %>%
    add_option(c("-i", "--input"),
               metavar = "DIRECTORY",
               help = "Path to raw MALDI results [required]") %>%
    add_option(c("-o", "--output"),
               metavar = "DIRECTORY",
               help = "Path to output directory [required]") %>%
    add_option("--min-mass",
               metavar = "NUMBER",
               default = 0,
               help = "Minimum mass to include in analysis [0]") %>%
    add_option("--max-mass",
               metavar = "NUMBER",
               default = Inf,
               help = "Maximum mass to include in analysis [Inf]")

  args <- parse_args2(parser)
  
  is_input_null <- is.null(args$options$input)
  is_output_null <- is.null(args$options$output)
  
  if (is_input_null || is_output_null) {
    print_help(parser)
    
    if (is_input_null) {
      cat("Required argument `--input` was not provided.\n")
    }
    if (is_output_null) {
      cat("Required argument `--output` was not provided.\n")
    }
    quit(save = "no", status = 1)
  }
}

main <- function() {
  args <- arguments()
  print(args)
}

load_spectra <- function(results_path, min_mass, max_mass) {
  spectra <-
    MALDIquantForeign::import(results_path, massRange = c(min_mass, max_mass))
  
  non_empty <-
    spectra %>%
    map_lgl(isEmpty) %>%
    all() %>%
    not()
  
  regular <-
    spectra %>%
    map_lgl(isRegular)
  
  if ((!non_empty) && regular) {
    cat("Input data are malformed\n")
    quit(save = "no", status = 1)
  }
  
  spectra
}

preprocess_spectra <- function(spectra) {
  transformed_spectra <-
    spectra %>%
    transformIntensity(method = "sqrt") %>%
    smoothIntensity(method = "SavitzkyGolay", halfWindowSize = 10)
  
  baseline_estimates <-
    map(transformed_spectra, function(spectrum) {
      estimateBaseline(spectrum, method = "SNIP", iterations = 100)
    })
  
  normalized_spectra <-
    transformed_spectra %>%
    removeBaseline(method = "SNIP", iterations = 100) %>%
    calibrateIntensity(method = "TIC") %>%
    alignSpectra(
      halfWindowSize = 20,
      SNR = 3,
      tolerance = 0.002,
      warpingMethod = "lowess"
    )
  
  normalized_spectra
}

average_technical_replicates <- function(spectra) {
  
  samples <-
    factor(map_chr(spectra, function(x) {
      metaData(x)$sampleName
    }))
  
  avg_spectra <-
    averageMassSpectra(spectra, labels = samples, method = "mean")
  
  avg_spectra
}

peak_detection <- function(spectrum, signal_to_noise_ratio) {
  
  noise <- estimateNoise(spectrum)
  
  peaks <-
    spectrum %>%
    detectPeaks(method = "MAD",
                halfWindowSize = 20,
                SNR = signal_to_noise_ratio) %>%
    binPeaks(tolerance = 0.002) %>%
    filterPeaks(minFrequency = 0.25)
  
  list("noise" = noise, "peaks" = peaks)
}

write_output <- function(spectra, peaks, noise) {
  
  
}
main()