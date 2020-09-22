# Massageplot[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

## Usage
![sashimiplot](example/sashimiplot.png)
```
Usage: sashimiplot [OPTIONS] COMMAND [ARGS]...

  A pure python sashimiplot, support the protein domain. Version: 1.0.0

Options:
  --help  Show this message and exit.

Commands:
  gene  Normal mode to generate sashimi plot
  junc  Junction mode, not need network to plot
  site  site mode, plot the last site coverage of the gene direction

```

## subcommand

### gene

add single gene model into your sashimiplot.

```
Usage: sashimiplot [OPTIONS] COMMAND [ARGS]...

  A pure python sashimiplot, support the protein domain. Version: 1.0.0

Options:
  --help  Show this message and exit.

Commands:
  gene  Normal mode to generate sashimi plot
  junc  Junction mode, not need network to plot
  site  site mode, plot the last site coverage of the gene direction

``` 

### junc

just for plot the sashimiplot

```
Usage: sashimiplot junc [OPTIONS]

  Junction mode, not need network to plot

Options:
  --gtf TEXT         The gtf file.
  --bam TEXT         Bam config file. There were two columns, label and file
                     path

  --fileout TEXT     The output name.
  --junc TEXT        The junction, it looks like chr:s:e
  --sj INTEGER       Only values greater than a certain value are displayed.
                     default: 1

  --pa TEXT          The pA site, if there were multiple sites, pls seperate
                     by `,`. default: None

  --wt TEXT          The weight for every pa sites from `--pa`, the number of
                     wt must be same to `--pa`,seperate by `,`

  --pa2 TEXT         The pA site with different color, if there were multiple
                     sites, pls seperate by `,`. default: None

  --wt2 TEXT         The weight for every pa sites from `--pa2`, the number of
                     wt must be same to `--pa`,seperate by `,`

  --focus TEXT       Highlight the given region. for one region: start-end, if
                     multiple, pls seperate by ,

  --ps TEXT          Library type, FR or RF.
  --peakfilter TEXT  peaks filter, default: None.
  --log TEXT         plot the log-tranformed expression values, and support
                     log2 and log10. default: None

  --prob TEXT        probability, given a region
  --id_keep TEXT     keep the given isoform, if there were multiple isoform
                     id, please seperate by comma

  --ssm TEXT         single-strand mode(ssm), given R1 or R2
  --ade              add expression for plot.
  --domain           add gene structure.
  --model TEXT       The deep learning model.
  --hl TEXT          highlight the given splice junction. Egg, sj1:sj2,sj3:sj4
  --inc TEXT         Only plot the given splice junction. Egg, sj1:sj2,sj3:sj4
  --dim TEXT         The picture's size,(width, height), default: 8,12
  --scale            Scale the count into 10%
  --bc TEXT          Barcode file for split single cell RNA seq data sets.
                     default: None

  --tag TEXT         Cell barcode and umi molecular tag in bam file, default:
                     CB,UB

  --help             Show this message and exit.

```

### site

plot the reads coverage based the most 3' site

```
Usage: sashimiplot site [OPTIONS]

  site mode, plot the last site coverage of the gene direction

Options:
  --gtf TEXT         The gtf file.
  --bam TEXT         Bam config file. There were two columns, label and file
                     path
  --fileout TEXT     The output name.
  --loc TEXT         The junction, it looks like chr:s:e
  --peakfilter TEXT  peaks filter, default: None.
  --verbose          set the logging level, if Ture -> INFO
  --log TEXT         plot the log-tranformed expression values, and support
                     log2 and log10. default: None
  --help             Show this message and exit.

```
