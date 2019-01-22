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
def gene(gtf, gene, bam, pa, fileout, offset, sj):
    if not all([gtf, gene, bam, fileout]):
        cli(['gene', '--help'])
        sys.exit(1)
    logging = set_logging("gene")
    wide = 8
    height = 12

    # geneinfo = Myinfo('ensembl.gene:{}'.format(gene),
    #                   'all',
    #                   'gene').loc
    #
    # mRNAobject = mRNA(geneinfo.chr,
    #                   geneinfo.start - offset,
    #                   geneinfo.end + offset,
    #                   gtf,
    #                   genename=gene
    #                   )
    '''
    1.21 add transcript support
    '''
    loggings.debug("prepare the mRNA data")
    mRNAobject = mRNA.gene(
        gene,
        gtf,
        offset
    )
    loggings.debug("retrieve expression data")
    bamdict = readbamlist(bam)
    bamlst = []
    for label, filepath in bamdict.items():
        bamlst.append({label: ReadDepth.determine_depth(filepath,
                                                        mRNAobject.chr,
                                                        mRNAobject.tstart,
                                                        mRNAobject.tend)})

    plot_density(bamlst,
                 mRNAobject,
                 fileout,
                 sj,
                 wide,
                 height,
                 pa
                 )

    # try:
    #     plot_density(bamlst,
    #                  mRNAobject,
    #                  fileout,
    #                  sj,
    #                  wide,
    #                  height,
    #                  pa
    #                  )
    # except:
    #     '''
    #     TODO need to add more information
    #     '''
    #     with open('failed.log', 'w') as faillog:
    #         if pa:
    #             faillog.write('{},{}'.format(gene, pa) + '\n')
    #         else:
    #             faillog.write(gene + '\n')


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
    loggings = set_logging("junc")

    if not all([gtf, bam, fileout, junc]):
        cli(['junc', '--help'])
        sys.exit(1)

    chr, s, e = junc.split(':')
    wide = 8
    height = 12
    pa = None
    domain = False
    # chr, tstart, tend, gtf
    loggings.debug("prepare the mRNA data")
    mRNAobject = mRNA(
        chr,
        s,
        e,
        gtf,
        exonstat=True
    )
    print(mRNAobject.tstart, mRNAobject.tend)
    bamdict = readbamlist(bam)
    bamlst = []
    loggings.debug("retrieve expression data")
    for label, filepath in bamdict.items():
        bamlst.append({label: ReadDepth.determine_depth(filepath,
                                                        mRNAobject.chr,
                                                        mRNAobject.tstart,
                                                        mRNAobject.tend)})
    loggings.debug("plot")
    plot_density(bamlst,
                 mRNAobject,
                 fileout,
                 sj,
                 wide,
                 height,
                 pa,
                 domain
                 )

    # try:
    #     plot_density(bamlst,
    #                  mRNAobject,
    #                  fileout,
    #                  sj,
    #                  wide,
    #                  height,
    #                  pa,
    #                  domain
    #                  )
    # except:
    #     '''
    #     TODO need to add more information
    #     '''
    #     with open('failed.log', 'w') as faillog:
    #         if pa:
    #             faillog.write('{},{}'.format(gene, pa) + '\n')
    #         else:
    #             faillog.write(gene + '\n')


cli.add_command(gene)
cli.add_command(junc)

if __name__ == '__main__':
    cli()
