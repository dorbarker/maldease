# Maldease

Maldease is a tool for rapidly processing [MALDI](https://en.wikipedia.org/wiki/Matrix-assisted_laser_desorption/ionization) mass spectra. It is currently in-development, but is usable in its current form. Maldease uses the excellent [MALDIquant](https://github.com/sgibb/maldiquant) library.

# Basic Usage

``` sh
maldease -i <path to directory of MALDI data> -o <output directory>
```

One or more defintions files may be placed in `data/`. Maldease will use these files for calling the presence or absence of peaks. The definition files are tab-delimited.

    target  mass    tolerance
    BoNT-A  2406    0.5
    BoNT-A  1203    1
    BoNT-A  1427    0.5
    BoNT-A  999 0.6

``` sh
$ maldease --help

Options:
    -h, --help
        Show this help message and exit

    -i DIRECTORY, --input=DIRECTORY
        Path to raw MALDI results [required]

    -o DIRECTORY, --output=DIRECTORY
        Path to output directory [required]

    --min-mass=NUMBER
        Minimum mass to include in analysis [0]

    --max-mass=NUMBER
        Maximum mass to include in analysis [Inf]

    --half-window-size=INT
        Half-window size for smoothing interval [20]

    --include-only=RANGES
        Consider only these mass ranges for peak reporting. 
              Lower and upper bounds are inclusive and separated by a dash. 
              Multiple ranges can be specified and are commas-delimited. 
              e.g. `--include-only '123-456,1000-2000'`

    --definitions=FILE
        Path to additional definitions for peak calling

    -v, --version
        Print the version and exit
```

# Data Citation

Example BoNT-A peaks are described in:

Kalb, *et al*. 2015. "Detection of Botulinum Toxins A, B, E, and F in Foods by Endopep-MS" <https://pubs.acs.org/doi/10.1021/jf505482b>
