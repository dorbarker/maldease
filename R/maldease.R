#!/usr/bin/env Rscript

for (pkg in c(
  "MALDIquant",
  "MALDIquantForeign",
  "optparse",
  "magrittr",
  "readr",
  "purrr",
  "here",
  "stringr",
  "tibble",
  "tidyr",
  "dplyr",
  "ggplot2",
  "writexl"
)) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

script_loc <-
    commandArgs(trailingOnly=FALSE) %>%
    keep(~str_detect(.x, "--file")) %>%
    str_remove("--file=") %>%
    tools::file_path_as_absolute()

project_loc <-
    script_loc %>%
    str_split(.Platform$file.sep) %>%
    unlist() %>%
    `[`(., 1:(length(.)-2)) %>%
    paste(collapse=.Platform$file.sep)

setwd(project_loc)

suppressMessages(here::i_am("R/maldease.R"))

arguments <- function() {

  parse_range <- function(range_str) {

    parsed_ranges <-
      range_str %>%
      str_split("\\s*,\\s*", simplify = TRUE) %>%
      str_split("-") %>%
      map(as.double) %>%
      map(sort)

    parsed_ranges
  }

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
               help = "Maximum mass to include in analysis [Inf]") %>%
    add_option(
      "--half-window-size",
      type = "integer",
      default = 20,
      metavar = "INT",
      help = "Half-window size for smoothing interval [20]"
    )  %>%
    add_option(
      "--include-only",
      default = "0-Inf",
      dest = "include_only_str",
      type = "character",
      metavar = "RANGES",
      help = "Consider only these mass ranges for peak reporting. \
              Lower and upper bounds are inclusive and separated by a dash. \
              Multiple ranges can be specified and are commas-delimited. \
              e.g. `--include-only '123-456,1000-2000'`"
    ) %>%
    add_option("--definitions",
               type = "character",
               metavar = "FILE",
               help = "Path to additional definitions for peak calling") %>%
    add_option(c("-v", "--version"),
               action = "store_true",
               default = FALSE,
               help = "Print the version and exit")


  args <- parse_args2(parser)

  if (args$option$version) {
    cat(paste("maldease", get_version()), "\n")
    quit(save = "no", status = 0)
  }

  args$options$include_only <- parse_range(args$options$include_only_str)

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

  analysis_time <- Sys.time()
  analysis_params <- format_analysis_parameters(args, analysis_time)

  spectra <-
    load_spectra(args$input, args$min_mass, args$max_mass, FALSE)

  target_definitions <- load_definitions(args$definitions)

  normalized_spectra <-
    preprocess_spectra(spectra, args$half_window_size)
  samples_avgSpectra <-
    average_technical_replicates(normalized_spectra)

  samples <- samples_avgSpectra$sample_names
  avg_spectra <- samples_avgSpectra$avg_spectra

  peaks_noise <-
    peak_detection_spectra(avg_spectra, 3, args$half_window_size)

  intensity_table <- get_intensity_table(avg_spectra,
                                         peaks_noise$peaks,
                                         samples,
                                         args$include_only)

  spectrum_plot <- draw_plots(avg_spectra, peaks_noise$peaks, args$include_only, analysis_time)

  calls <- call_positives(intensity_table, target_definitions)

  write_output(spectrum_plot, intensity_table, analysis_params, calls, args$output)

}

get_version <- function() {
  version <-
    read_lines(here("DESCRIPTION")) %>%
    keep(~ grepl("Version", .x)) %>%
    str_split(": ") %>%
    unlist() %>%
    .[2]

  version
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

load_definitions <- function(additional_path) {

  # definition file format (no headers, comma-delimited):
  # name1,mass,mass_tolerance
  # name1,mass,mass_tolerance
  # name2,mass,mass_tolerance
  # etc

  definitions_files <-
    here("data") %>%
    list.files(full.names = TRUE)

  if (!is.null(additional_path)) {
    definitions_files <- c(definitions_files, additional_path)
  }

  definitions <-
    definitions_files %>% map_dfr(~ read_tsv(.x, col_types = "cdd"))

  definitions
}

preprocess_spectra <- function(spectra, half_window_size) {
  transformed_spectra <-
    spectra %>%
    # transformIntensity(method = "sqrt") %>%
    suppressWarnings({
      smoothIntensity(., method = "SavitzkyGolay",
                      halfWindowSize = half_window_size)
    })


  baseline_estimates <-
    map(transformed_spectra, function(spectrum) {
      estimateBaseline(spectrum, method = "SNIP", iterations = 100)
    })
  normalized_spectra <-
    transformed_spectra %>%
    removeBaseline(method = "SNIP", iterations = 100) %>%
    # calibrateIntensity(method = "TIC") %>%
    alignSpectra(
      halfWindowSize = half_window_size,
      SNR = 3,
      tolerance = 0.02,
      noiseMethod="SuperSmoother",
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

  list(
    "sample_names" = map_chr(avg_spectra, ~ metaData(.x)$sampleName),
    "avg_spectra" = avg_spectra
  )
}

peak_detection_spectrum <-  function(spectrum,
                                     signal_to_noise_ratio,
                                     half_window_size) {
    peaks <-
      spectrum %>%
      detectPeaks(method = "SuperSmoother",
                  halfWindowSize = half_window_size,
                  SNR = signal_to_noise_ratio) %>%
      binPeaks(method="relaxed", tolerance = 0.002) %>%
      filterPeaks(minFrequency = 0.25)

    peaks
  }

peak_detection_spectra <- function(spectra,
                                   signal_to_noise_ratio,
                                   half_window_size) {
    noises <- map(spectra, ~ estimateNoise(.x, method = "SuperSmoother"))

    peaks <-
      detectPeaks(spectra,
                  noiseMethod = "SuperSmoother",
                  halfWindowSize = half_window_size,
                  SNR = signal_to_noise_ratio) %>%
      binPeaks(tolerance = 0.002) %>%
      filterPeaks(minFrequency = 0.25)

    list("noise" = noises, "peaks" = peaks)
  }

filter_mass_on_ranges <- function(mass_spec_table, incl_ranges) {

  map_dfr(incl_ranges, ~ filter(mass_spec_table, between(mass, .x[1], .x[2])))
}

get_intensity_table <- function(spectra, peaks, samples, incl_ranges) {
  m <- intensityMatrix(peaks, spectra)
  mz <- attr(m, "mass")

  mass_snr <-
    map_dfr(peaks,
            ~ tibble(
              sample = .x@metaData$sampleName,
              mass = .x@mass,
              snr = .x@snr
            ))

  intensity_table <-
    t(m) %>%
    as_tibble(.name_repair = "minimal") %>%
    set_colnames(samples) %>%
    mutate(mass = mz) %>%
    pivot_longer(-mass,
                 names_to = "sample",
                 values_to = "intensity") %>%
    left_join(mass_snr, by = c("mass", "sample")) %>%
    mutate(noise = intensity / snr) %>%
    drop_na() %>%
    filter_mass_on_ranges(incl_ranges)

  intensity_table
}

draw_plots <- function(spectra, peaks, incl_ranges, analysis_time) {

  peak_tables <-
    map_dfr(
      peaks,
      ~ tibble(
        sample = .x@metaData$sampleName,
        mass = .x@mass,
        peak_intensity = .x@intensity
      )
    ) %>%
    filter_mass_on_ranges(incl_ranges)

  spectrum_tables <- imap_dfr(spectra,
                              ~ tibble(
                                sample = .y,
                                mass = .x@mass,
                                intensity = .x@intensity)
                              )

  p <-
    ggplot(spectrum_tables, aes(x = mass, y = intensity)) +
    geom_line() +
    geom_point(
      data = peak_tables,
      aes(x = mass, y = peak_intensity),
      shape = 4,
      colour = "red"
    ) +
    labs(title = analysis_time) +
    facet_wrap( ~ sample,
                ncol = 1,
                scales = "fixed") +
    theme_linedraw()

  p
}

format_analysis_parameters <- function(arg, analysis_time) {

  analysis_time_str <- paste("analysis_time", analysis_time, sep = ": ")
  cli_args <- paste(names(arg), arg, sep = ": ", collapse = "\n")

  analysis_parameters_log <- paste(analysis_time_str, cli_args, "", sep = "\n")

  analysis_parameters_log
}

call_positives <- function(intensity_table, target_definitions) {

  target_definitions %>%
    mutate(
      lower = mass - tolerance,
      upper = mass + tolerance,
    ) %>%
    rowwise() %>%
    mutate(
      positive = between(intensity_table$mass, lower, upper) %>% any()
    ) %>%
    select(target, mass, positive)
}

write_output <- function(spectrum_plot, intensity_table, params, calls, outpath) {
  if (!dir.exists(outpath)) {
    dir.create(outpath)
  }

  plot_path <- file.path(outpath, "mass_spectrum.png")
  table_path <- file.path(outpath, "intensity_table.tsv")
  excel_path <- file.path(outpath, "intensity_table.xlsx")
  params_path <- file.path(outpath, "params.txt")
  calls_path <- file.path(outpath, "calls.tsv")

  ggsave(
    filename = plot_path,
    plot = spectrum_plot,
    height = 15,
    width = 30,
    units = "cm"
  )

  write_tsv(intensity_table, table_path)
  write_tsv(calls, calls_path)
  write_xlsx(intensity_table, excel_path)
  write_lines(params, params_path)
}

main()
