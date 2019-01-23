# Massageplot[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

## Usage
![sashimiplot](example/sashimiplot.png)
```
Usage: sashimiplot [OPTIONS] COMMAND [ARGS]...

  A pure python sashimiplot, support the protein domain. Version: 1.1.0

Options:
  --help  Show this message and exit.

Commands:
  gene
  junc
```

## subcommand

### gene
add single gene model into your sashimiplot.

```
Usage: sashimiplot gene [OPTIONS]

Options:
  --gtf TEXT      The gtf file.
  --gene TEXT     The ensembl gene id only for now.
  --bam TEXT      Bam config file. There were two columns, label and file path
  --pa INTEGER    The pA site.
  --fileout TEXT  The output name.
  --help          Show this message and exit.
``` 

### junc

just for plot the sashimiplot
```
Usage: sashimiplot junc [OPTIONS]

  Junction mode, not need network to plot

Options:
  --gtf TEXT      The gtf file.
  --bam TEXT      Bam config file. There were two columns, label and file path
  --fileout TEXT  The output name.
  --junc TEXT     The junction, it looks like chr:s:e
  --sj INTEGER    Only values greater than a certain value are displayed.
                  default: 1
  --help          Show this message and exit.
```
