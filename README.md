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

  Normal mode to generate sashimi plot

Options:
  --gtf TEXT        The gtf file.
  --gene TEXT       The ensembl gene id only for now.
  --bam TEXT        Bam config file. There were two columns, label and file
                    path
  --pa INTEGER      The pA site. default: None
  --fileout TEXT    The output name.
  --offset INTEGER  Add an offset number to broaden the interval. default: 0
  --sj INTEGER      The min splice jucntion count to show. default: 1
  --focus TEXT      Highlight the given region. start-end
  --help            Show this message and exit.

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
  --pa TEXT       The pA site, if there were multiple sites, pls seperate by
                  `,`. default: None
  --focus TEXT    Highlight the given region. for one region: start-end, if
                  multiple, pls seperate by ,
  --ps            plot the site coverage, default: False.
  --help          Show this message and exit.
```

### site

plot the reads coverage based the most 3' site

```
Usage: sashimiplot site [OPTIONS]

  site mode, plot the last site coverage of the gene direction

Options:
  --gtf TEXT      The gtf file.
  --bam TEXT      Bam config file. There were two columns, label and file path
  --fileout TEXT  The output name.
  --loc TEXT      The junction, it looks like chr:s:e
  --help          Show this message and exit.
```
