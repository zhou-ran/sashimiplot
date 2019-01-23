#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/1/16 4:07 PM
__author__ = 'Zhou Ran'
import sys
import click

from .helper import set_logging
from .ReadDepth import ReadDepth
from .mRNA import mRNA
from .plot import plot_density
from .utils import readbamlist

logger = set_logging("MAIN")


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
              help="The pA site. default: None"
              )
@click.option('--fileout',
              type=str,
              help="The output name."
              )
@click.option('--offset',
              type=int,
              default=0,
              help='Add an offset number to broaden the interval. default: 0')
@click.option('--sj',
              type=int,
              default=1,
              help='The min splice jucntion count to show. default: 1')
@click.option('--focus',
              default=None,
              help="Highlight the given region. start-end")
def gene(gtf, gene, bam, pa, fileout, offset, sj, focus):
    if not all([gtf, gene, bam, fileout]):
        cli(['gene', '--help'])
        sys.exit(1)

    wide = 8
    height = 12

    if focus:
        focus = focus.split('-')

    '''
    1.21 add transcript support
    '''
    logger.info("prepare the mRNA data")
    mRNAobject = mRNA.gene(
        gene,
        gtf,
        offset
    )

    logger.info("retrieve expression data")
    bamdict = readbamlist(bam)
    bamlst = []
    for label, filepath in bamdict.items():
        bamlst.append({label: ReadDepth.determine_depth(filepath,
                                                        mRNAobject.chr,
                                                        mRNAobject.tstart,
                                                        mRNAobject.tend)})
    logger.info("plot")
    try:
        plot_density(bamlst,
                     mRNAobject,
                     fileout,
                     sj,
                     wide,
                     height,
                     pa
                     )

    except Exception as e:
        logger.error(e)


@click.command()
@click.option('--gtf',
              type=str,
              help="The gtf file."
              )
@click.option('--bam',
              type=str,
              help="Bam config file. There were two columns, label and file path"
              )
@click.option('--fileout',
              type=str,
              help="The output name."
              )
@click.option('--junc',
              help="The junction, it looks like chr:s:e"
              )
@click.option('--sj',
              type=int,
              default=1,
              help="Only values greater than a certain value are displayed. default: 1")
def junc(gtf, bam, fileout, junc, sj):
    """
    Junction mode, not need network to plot
    :return:
    """

    if not all([gtf, bam, fileout, junc]):
        cli(['junc', '--help'])
        sys.exit(1)

    chr, s, e = junc.split(':')
    wide = 8
    height = 12
    pa = None
    domain = False
    # chr, tstart, tend, gtf
    logger.info("prepare the mRNA data")
    mRNAobject = mRNA(
        chr,
        s,
        e,
        gtf,
        exonstat=True
    )

    bamdict = readbamlist(bam)
    bamlst = []
    logger.info("retrieve expression data")
    for label, filepath in bamdict.items():
        bamlst.append({label: ReadDepth.determine_depth(filepath,
                                                        mRNAobject.chr,
                                                        mRNAobject.tstart,
                                                        mRNAobject.tend)})
    logger.info("plot")
    try:
        plot_density(bamlst,
                     mRNAobject,
                     fileout,
                     sj,
                     wide,
                     height,
                     pa,
                     domain
                     )
    except Exception as e:
        logger.error(e)


cli.add_command(gene)
cli.add_command(junc)

if __name__ == '__main__':
    cli()
