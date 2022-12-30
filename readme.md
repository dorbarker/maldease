# Maldease

Maldease is a tool for rapidly processing [MALDI](https://en.wikipedia.org/wiki/Matrix-assisted_laser_desorption/ionization) mass spectra. It is currently in-development, but is usable in its current form. Maldease uses the excellent [MALDIquant](https://github.com/sgibb/maldiquant) library.

# Basic Usage

``` sh
maldease --negative-control=<path> --output=<path> --definitions=<tsv> <experiments>...
```

```
Usage:
  maldease --negative-control=<path> --output=<path> --definitions=<tsv>... [--peaks-only|--include-only=<ranges>] [options] <experiments>...
  maldease --generate-example-definitions
  maldease -v | --version
  maldease -h | --help

Options:
  -h --help                            Print this help and exit
  -v --version                         Print program version and exit
  -n <path> --negative-control=<path>  Path to negative control MALDI results
  -o <path> --output=<path>            Output directory path
  --half-window-size=<int>             Half-window size for smoothing and peak
                                       detection [default: 20]
  --include-only=<ranges>              Consider only these mass ranges for
                                       peak reporting. Lower and upper bounds
                                       are inclusive and separated by a dash.
                                       Multiple ranges can be specified and
                                       are comma delimited,
                                       e.g.--include-only '123-456,1000-2000'
                                       [default: 0-Inf]
  --peaks-only                         Only report peaks within the tolerance
                                       specified by --definitions
  -d <tsv> --definitions=<tsv>         One or more paths to tab-separated
                                       defintions of target peaks
  --generate-example-definitions       Generate an example definitions file
                                       and write it to stdout 
```

# Data Citation

Example BoNT-A peaks are described in:

Kalb, *et al*. 2015. "Detection of Botulinum Toxins A, B, E, and F in Foods by Endopep-MS" <https://pubs.acs.org/doi/10.1021/jf505482b>
