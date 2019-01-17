#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/1/16 4:07 PM
__author__ = 'Zhou Ran'
# import sys
import click

from .ReadDepth import ReadDepth
from .FetchGene import Myinfo
from .mRNA import mRNA
from .plot import plot_density
from .utils import readbamlist


@click.group()
def cli():
    """A pure python sashimiplot, support the protein domain. Version: 1.0.0"""
    pass


@click.command()
@click.option('--gtf',
              type=str,
              help="The gtf file."
              )
@click.option('--gene',
              type=str,
              help="The ensembl gene id only for now."
              )
@click.option('--bam',
              type=str,
              help="Bam config file. There were two columns, label and file path"
              )
@click.option('--pa',
              type=int,
              default=None,
              help="The pA site."
              )
@click.option('--fileout',
              type=str,
              help="The output name."
              )
def gene(gtf, gene, bam, pa, fileout):
    click.echo("Plot a sashimi in the given gene!")
    geneinfo = Myinfo('ensembl.gene:{}'.format(gene),
                      'all',
                      'gene').loc

    mRNAobject = mRNA(geneinfo.chr,
                      geneinfo.start,
                      geneinfo.end,
                      gtf,
                      genename=gene)

    bamdict = readbamlist(bam)
    bamlst = []
    for label, filepath in bamdict.items():
        bamlst.append({label: ReadDepth.determine_depth(filepath,
                                                        geneinfo.chr,
                                                        geneinfo.start,
                                                        geneinfo.end)})

    plot_density(bamlst,
                 mRNAobject,
                 "+" if geneinfo.strand > 0 else "-",
                 fileout,
                 pa)


@click.command()
def junc():
    click.echo("Plot a sashimi in the given interval")


cli.add_command(gene)
cli.add_command(junc)

if __name__ == '__main__':
    cli()
