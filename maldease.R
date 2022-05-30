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
              "purrr",
              "tibble",
              "tidyr",
              "dplyr",
              "ggplot2")) {
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
  args
}

main <- function() {
  args <- arguments()$options

  spectra <- load_spectra(args$input, args$min_mass, args$max_mass, FALSE)
  normalized_spectra <- preprocess_spectra(spectra)
  samples_avgSpectra <- average_technical_replicates(spectra)
  
  samples <- samples_avgSpectra$sample_names
  avg_spectra <- samples_avgSpectra$avg_spectra
  
  peaks_noise <- peak_detection_spectra(avg_spectra, 3)
  
  intensity_table <- get_intensity_table(peaks_noise$peaks, peaks_noise$noise, avg_spectra, samples)
  
  intensity_table %>% format_tsv() %>% cat()
}

load_spectra <- function(results_path, min_mass, max_mass, verbose) {
  
  spectra <-
    if (verbose) {
      MALDIquantForeign::import(results_path,
                                massRange = c(min_mass, max_mass))
    } else {
      suppressWarnings({
        MALDIquantForeign::import(results_path,
                                  massRange = c(min_mass, max_mass))
      })
    }

  
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
    #calibrateIntensity(method = "TIC") %>%
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
  
  list("sample_names" = map_chr(avg_spectra, ~ metaData(.x)$sampleName),
       "avg_spectra" = avg_spectra)
}

peak_detection_spectrum <- function(spectrum, signal_to_noise_ratio) {
  
  peaks <-
    spectrum %>%
    detectPeaks(method = "MAD",
                halfWindowSize = 20,
                SNR = signal_to_noise_ratio) %>%
    binPeaks(tolerance = 0.002) %>%
    filterPeaks(minFrequency = 0.25)
  
  peaks
}

peak_detection_spectra <- function(spectra, signal_to_noise_ratio) {
  
  noises <- map(spectra, estimateNoise)
  #peaks <- map(spectra, ~ peak_detection_spectrum(.x, signal_to_noise_ratio))
  peaks <-
    detectPeaks(spectra,
                method = "MAD",
                halfWindowSize = 20,
                SNR = signal_to_noise_ratio) %>% 
    binPeaks(tolerance = 0.002) %>% 
    filterPeaks(minFrequency = 0.25)
  
  
  list("noise" = noises, "peaks" = peaks)
}

get_intensity_table <- function(peaks, noise, spectra, samples) {
  
  m <- intensityMatrix(peaks, spectra)
  mz <- attr(m, "mass")
  
  intensity_table <- 
    t(m) %>% 
    as_tibble() %>% 
    set_colnames(samples) %>% 
    mutate(mass = mz) %>% 
    pivot_longer(-mass, names_to = "sample", values_to = "intensity")
    
  noises <- tibble()
  for (n in names(noise)) {
    cur_noise <-
      noise[[n]] %>%
      as_tibble() %>%
      # rename(noise = intensity) %>%
      transmute(sample = n, noise = intensity) %>%
      sample_n(1)
      

    noises <- bind_rows(noises, cur_noise)

  }
  
  intensity_table_with_snr <-
    intensity_table %>% 
    left_join(noises) %>% 
    mutate(snr = intensity / noise)
}

write_output <- function(spectra, peaks, noise) {

  for (i in seq_along(spectra)) {
    
  }  
  
}
main()