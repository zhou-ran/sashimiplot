# Massageplot

## Usage
```
Usage: sashimiplot [OPTIONS] COMMAND [ARGS]...

  A pure python sashimiplot, support the protein domain. Version: 1.0.0

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
