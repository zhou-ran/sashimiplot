#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/1/16 4:07 PM

__author__ = 'Zhou Ran'
import sys

import click
from loguru import logger

from .ReadDepth import ReadDepth
from .mRNA import mRNA
from .plot import plot_density
from .sitedepth import SiteDepth
from .siteplot import plot_density_site
from .utils import readbamlist

logger.remove(0)
logger.add('sashimiplot.log', rotation='10 MB', colorize=True,
           level="DEBUG")

logger.add(sys.stderr, level="INFO", colorize=True)


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
@click.option('--verbose',
              is_flag=True,
              help='set the logging level, if Ture -> INFO')
@click.option('--log',
              default=None,
              type=str,
              help="plot the log-tranformed expression values, and support log2 and log10. default: None"
              )
def gene(gtf, gene, bam, pa, fileout, offset, sj, focus, log, verbose):
    """
    Normal mode to generate sashimi plot
    """

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
    bamdict, colordict = readbamlist(bam)
    bamlst = []
    # colors = []
    for label, filepath in bamdict.items():

        bamlst.append({label: ReadDepth.determine_depth(filepath,
                                                        mRNAobject.chr,
                                                        mRNAobject.tstart,
                                                        mRNAobject.tend)})
        # colors.append(color)

    logger.info("plot")
    try:
        plot_density(bamlst,
                     mRNAobject,
                     fileout,
                     sj,
                     wide,
                     height,
                     colors = colordict,
                     pasite=pa,
                     focus=focus,
                     logtrans=log
                     )

    except Exception as e:
        logger.error("Error information found in {}, pls check the gene command".format(gene))
        logger.exception(e)


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
@click.option('--pa',
              type=str,
              default=None,
              help="The pA site, if there were multiple sites, pls seperate by `,`. default: None"
              )
@click.option('--wt',
              type=str,
              default=None,
              help="The weight for every pa sites from `--pa`, the number of wt must be same to `--pa`,seperate by `,`"
              )

@click.option('--pa2',
              type=str,
              default=None,
              help="The pA site with different color, if there were multiple sites, pls seperate by `,`. default: None"
              )
@click.option('--wt2',
              type=str,
              default=None,
              help="The weight for every pa sites from `--pa2`, the number of wt must be same to `--pa`,seperate by `,`"
              )

@click.option('--focus',
              default=None,
              help="Highlight the given region. for one region: start-end, if multiple, pls seperate by ,")
@click.option('--ps',
              default=None,
              type=str,
              help="Library type, FR or RF."
              )
@click.option('--peakfilter',
              default=None,
              type=str,
              help="peaks filter, default: None."
              )
@click.option('--log',
              default=None,
              type=str,
              help="plot the log-tranformed expression values, and support log2 and log10. default: None"
              )
@click.option('--prob',
              default=None,
              type=str,
              help="probability, given a region"
              )
@click.option('--id_keep',
              default=None,
              type=str,
              help="keep the given isoform, if there were multiple isoform id, please seperate by comma"
              )

@click.option('--ssm',
              default=None,
              type=str,
              help="single-strand mode(ssm), given R1 or R2"
              )
@click.option('--ade',
              is_flag=True,
              help="add expression for plot."
              )
@click.option('--domain',
              is_flag=True,
              help="add gene structure."
              )
@click.option('--model',
              default=None,
              type=str,
              help="The deep learning model.")

@click.option('--hl',
              default=None,
              type=str,
              help="highlight the given splice junction.egg, sj1:sj2,sj3:sj4")
@click.option('--dim',
              default=None,
              type=str,
              help="The picture's size,(width, height), default: 8,12")
@click.option('--scale',
              is_flag=True,
              help="Scale the count into 10%")
# @click.option('--verbose',
#               is_flag=True,
#               help='set the logging level, if Ture -> INFO')
def junc(gtf,
         bam,
         fileout,
         junc,
         sj,
         pa,
         wt,
         pa2,
         wt2,
         focus,
         ps,
         peakfilter,
         log,
         prob,
         ssm,
         domain,
         hl,
         dim,
         # verbose,
         model,
         ade,
         id_keep,
         scale
         ):
    """
    Junction mode, not need network to plot
    """

    if not all([gtf, bam, fileout, junc]):
        cli(['junc', '--help'])
        sys.exit(1)

    chr, s, e = junc.split(':')
    # print(scale)
    logger.info("prepare the mRNA data")
    mRNAobject = mRNA(
        chr,
        s,
        e,
        gtf,
        exonstat=True if not domain else False
    )

    # focus the given regions
    if focus:
        focus = focus.split(',')
        focus = map(lambda x: x.split('-'), focus)

    # highlight the splice junction
    if hl:
        # incase somebody add space in the given splice junction
        hl = set(map(lambda x: x.strip(), hl.split(',')))

    # id to keep
    if id_keep:
        id_keep = set(id_keep.split(','))

    # figdim
    if dim:
        width, height = map(lambda x: float(x), map(lambda x: x.strip(), dim.split(',')))
    else:
        width, height = 8.0, 12.0

    bamdict, colordict = readbamlist(bam)
    bamlst = []
    bamsitelst = [] if ps else None
    logger.info("retrieve expression data")

    for label, filepath in bamdict.items():
        readdepth_ = ReadDepth.generateobj()

        for bam_ in filepath:
            readdepth_ += ReadDepth.determine_depth(bam_,
                                                    mRNAobject.chr,
                                                    mRNAobject.tstart,
                                                    mRNAobject.tend,
                                                    scale=scale,
                                                    readFilter=peakfilter)
        bamlst.append({label: readdepth_})
        if ps:
            sitedepth_ = SiteDepth.generateobj()
            for bam_ in filepath:
                sitedepth_ += SiteDepth.determine_depth(bam_,
                                                        mRNAobject.chr,
                                                        mRNAobject.tstart,
                                                        mRNAobject.tend,
                                                        ps,
                                                        singlestrand=ssm,
                                                        readFilter=peakfilter)
            bamsitelst.append({label: sitedepth_})
    logger.info("plot")

    try:
        plot_density(bamlst,
                     mRNAobject,
                     fileout,
                     sj,
                     width=width,
                     height=height,
                     colors = colordict,
                     pasite=pa,
                     wt_pasite=wt,
                     pasite2=pa2,
                     wt_pasite2=wt2,
                     focus=focus,
                     domain=domain,
                     sitedepth=bamsitelst,
                     logtrans=log,
                     prob=prob,
                     model=model,
                     addexpress=ade,
                     id_keep=id_keep,
                     hl=hl
                     )
    except Exception as e:
        logger.error("Error information found in {}, pls check the splicing region".format(junc))
        logger.exception(e)


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
@click.option('--loc',
              help="The junction, it looks like chr:s:e"
              )
@click.option('--peakfilter',
              default=None,
              type=str,
              help="peaks filter, default: None."
              )
@click.option('--verbose',
              is_flag=True,
              help='set the logging level, if Ture -> INFO')
@click.option('--log',
              default=None,
              type=str,
              help="plot the log-tranformed expression values, and support log2 and log10. default: None"
              )
def site(gtf,
         bam,
         fileout,
         loc,
         peakfilter,
         verbose,
         log
         ):
    """
    site mode, plot the last site coverage of the gene direction
    """

    if not all([gtf, bam, fileout, loc]):
        cli(['site', '--help'])
        sys.exit(1)

    chr, s, e = loc.split(':')

    logger.info("prepare the mRNA data")
    mRNAobject = mRNA(
        chr,
        s,
        e,
        gtf,
        exonstat=True
    )

    bamdict, colordict = readbamlist(bam)
    bamlst = []
    logger.info("retrieve expression data")

    for label, filepath in bamdict.items():
        readdepth_ = ''
        for bam_ in filepath:
            if readdepth_ == '':
                readdepth_ = SiteDepth.determine_depth(bam_,
                                                       mRNAobject.chr,
                                                       mRNAobject.tstart,
                                                       mRNAobject.tend,
                                                       "FR",
                                                       readFilter=peakfilter)
            else:
                readdepth_ += SiteDepth.determine_depth(bam_,
                                                        mRNAobject.chr,
                                                        mRNAobject.tstart,
                                                        mRNAobject.tend,
                                                        "FR",
                                                        readFilter=peakfilter)

        bamlst.append({label: readdepth_})

    logger.info("plot")
    try:
        plot_density_site(bamlst,
                          mRNAobject,
                          fileout=fileout,
                          logtrans=log
                          )
    except Exception as e:
        logger.error("Error information found in {}, pls check the splicing region".format(junc))
        logger.exception(e)

cli.add_command(gene)
cli.add_command(junc)
cli.add_command(site)

if __name__ == '__main__':
    cli()
