#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/1/16 4:07 PM
from functools import reduce
# from operator import add

__author__ = 'Zhou Ran'
import sys

import click
from loguru import logger
from collections import defaultdict
from random import sample as random_sample

from .ReadDepth import ReadDepth
from .mRNA import mRNA
from .plot import plot_density
from .sitedepth import SiteDepth
from .siteplot import plot_density_site
from .utils import read_track, readbamlist
from .Constant import COLOR
from .annot_track import AnnotTrack

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
@click.option('--dim',
              default=None,
              type=str,
              help="The picture's size,(width, height), default: 8,12")
@click.option('--inc',
              default=None,
              type=str,
              help="Only plot the given splice junction. Egg, sj1:sj2,sj3:sj4")
@click.option('--ie',
              default='15,1',
              type=str,
              help="Intron and exon scale for plot the mRNA track, default: 15,1")
@click.option('--log',
              default=None,
              type=str,
              help="plot the log-tranformed expression values, and support log2 and log10. default: None"
              )
def gene(gtf, gene, bam, pa, fileout, offset, sj, focus, log, dim, inc, ie, verbose):
    """
    Normal mode to generate sashimi plot
    """

    if not all([gtf, gene, bam, fileout]):
        cli(['gene', '--help'])
        sys.exit(1)

    # figdim
    if dim:
        width, height = map(lambda x: float(x), map(lambda x: x.strip(), dim.split(',')))
    else:
        width, height = 8.0, 12.0

    if focus:
        focus = focus.split('-')
    # which sj will be plot
    if inc:
        # incase somebody add space in the given splice junction
        inc = list(map(lambda x: x.strip(), inc.split(',')))
    try:
        intron_scale, exon_scale = map(int, map(lambda x: x.strip(), ie.split(',')))
    except ValueError as e:
        logger.error('Error intron and exon scale information.')
        logger.exception(e)

    '''
    1.21 add transcript support
    '''
    logger.info("prepare the mRNA data")
    mRNAobject = mRNA.gene(
        gene,
        gtf,
        offset
    )

    bamdict, colordict = readbamlist(bam)
    bamlst = []
    logger.info("retrieve expression data")

    for label, filepath in bamdict.items():
        readdepth_ = ReadDepth.generateobj()

        for bam_ in filepath:
            readdepth_ += ReadDepth.determine_depth(bam_,
                                                    mRNAobject.chr,
                                                    mRNAobject.tstart,
                                                    mRNAobject.tend)
        bamlst.append({label: readdepth_})

    logger.info("plot")
    try:
        plot_density(bamlst,
                     mRNAobject,
                     fileout,
                     sj,
                     width=width,
                     height=height,
                     colors=colordict,
                     pasite=pa,
                     focus=focus,
                     include_sj=inc,
                     logtrans=log,
                     intron_scale=intron_scale,
                     exon_scale=exon_scale
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
@click.option('--trackline',
              type=str,
              default=None,
              help="The track line on annotation, if there were multiple sites, pls seperate by `,`. default: None"
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
@click.option('--sameymax',
              is_flag=True,
              help="same or not same ymax."
              )
@click.option('--model',
              default=None,
              type=str,
              help="The deep learning model.")
@click.option('--hl',
              default=None,
              type=str,
              help="highlight the given splice junction. Egg, sj1:sj2,sj3:sj4")
@click.option('--inc',
              default=None,
              type=str,
              help="Only plot the given splice junction. Egg, sj1:sj2,sj3:sj4")
@click.option('--dim',
              default=None,
              type=str,
              help="The picture's size,(width, height), default: 8,12")
@click.option('--scale',
              is_flag=True,
              help="Scale the count into 10%")
@click.option('--bc',
              default=None,
              type=str,
              help="Barcode file for split single cell RNA seq data sets. default: None")
@click.option('--tag',
              default='CB,UB',
              type=str,
              help="Cell barcode and umi molecular tag in bam file, default: CB,UB")
@click.option('--ie',
              default='15,1',
              type=str,
              help="Intron and exon scale for plot the mRNA track, default: 15,1")
@click.option('--co',
              default=None,
              type=str,
              help="Cluster orders.")
@click.option('--track',
              default=None,
              type=str,
              help="Additional track file, like miRNA binding sites as bed formats.")
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
         inc,
         dim,
         model,
         ade,
         id_keep,
         scale,
         bc,
         ie,
         tag,
         co,
         trackline,
         track,
         sameymax
         ):
    """
    Junction mode, not need network to plot
    """

    if not all([gtf, bam, fileout, junc]):
        cli(['junc', '--help'])
        sys.exit(1)

    chrom, s, e = junc.split(':')

    logger.info("prepare the mRNA data")
    mRNAobject = mRNA(
        chrom,
        s,
        e,
        gtf,
        exonstat=True if not domain else False
    )
    if track:
        track = AnnotTrack(
            read_track(trackfile=track),
            mRNAobject.chr,
            mRNAobject.tstart,
            mRNAobject.tend
        ).tracklist

    if bc:
        cb_tag, umi_tag = map(lambda x: x.strip(), tag.split(','))
        cell_cluster = set()
        # sample_names = set()

        sample_cell_cluster = defaultdict(lambda: defaultdict())
        with open(bc) as bc_fh:
            for line in bc_fh:
                if not line:
                    continue
                cell_name, cluster_info = line.strip().split('\t')
                sample_name, cell_id = cell_name.split('_')
                cell_id = cell_id.split('-')[0]

                sample_cell_cluster[sample_name][cell_id] = cluster_info

                cell_cluster.add(cluster_info)
                # sample_names.add(sample_name)
    try:
        intron_scale, exon_scale = map(int, map(lambda x: x.strip(), ie.split(',')))
    except ValueError as e:
        logger.error('Error intron and exon scale information.')
        logger.exception(e)

    # focus the given regions
    if focus:
        focus = focus.split(',')
        focus = map(lambda x: x.split('-'), focus)

    # highlight the splice junction
    if hl:
        # incase somebody add space in the given splice junction
        hl = list(map(lambda x: x.strip(), hl.split(',')))

    # which sj will be plot
    if inc:
        # incase somebody add space in the given splice junction
        inc = list(map(lambda x: x.strip(), inc.split(',')))

    # id to keep
    if id_keep:
        id_keep = set(id_keep.split(','))

    # fig dim
    if dim:
        width, height = map(lambda x: float(x), map(lambda x: x.strip(), dim.split(',')))
    else:
        width, height = 8.0, 12.0

    bam_dic, color_dic = readbamlist(bam)
    bam_cov = {}
    bam_site_cov = {} if ps else None
    logger.info("retrieve expression data")
    tmp_list = []

    for label, filepath in bam_dic.items():
        if not bc:
            read_depth = ReadDepth.generateobj()
            for bam_ in filepath:
                read_depth += ReadDepth.determine_depth(
                    bam_,
                    mRNAobject.chr,
                    mRNAobject.tstart,
                    mRNAobject.tend,
                    scale=scale,
                    readFilter=peakfilter
                )

            bam_cov[label] = read_depth

            if ps:
                site_depth = SiteDepth.generateobj()
                for bam_ in filepath:
                    site_depth += SiteDepth.determine_depth(
                        bam_,
                        mRNAobject.chr,
                        mRNAobject.tstart,
                        mRNAobject.tend,
                        ps,
                        singlestrand=ssm,
                        readFilter=peakfilter
                    )
                bam_site_cov[label] = site_depth
        else:
            for bam_ in filepath:
                tmp_list.append(
                    ReadDepth.determine_depth(
                        bam_,
                        mRNAobject.chr,
                        mRNAobject.tstart,
                        mRNAobject.tend,
                        barcode_info=sample_cell_cluster[label],
                        cell_tag=cb_tag,
                        umi_tag=umi_tag,
                        scale=scale,
                        readFilter=peakfilter
                    )
                )

            if ps:
                # single cell bam, no site coverage. not available
                bam_site_cov = None
    if bc:
        for cluster in cell_cluster:
            cluster_use_obj = []
            for sample in tmp_list:

                try:
                    cluster_use_obj.append(
                        sample[cluster]
                    )
                except KeyError:
                    continue

            bam_cov[cluster] = reduce(ReadDepth.__add__, cluster_use_obj)

    if co:
        order_dic = {}
        try:
            with open(co) as co_fh:
                for line in co_fh:
                    if not line:
                        continue
                    try:
                        id_use, col_use = line.strip().split('\t')
                        order_dic[id_use] = col_use
                    except ValueError:
                        order_dic[id_use] = COLOR[random_sample(range(len(COLOR)), 1)]
        except Exception as e:
            logger.exception(e)
    else:
        order_dic = None

    logger.info("plot")

    try:
        plot_density(bam_cov,
                     mRNAobject,
                     fileout,
                     sj,
                     width=width,
                     height=height,
                     colors=color_dic,
                     pasite=pa,
                     wt_pasite=wt,
                     pasite2=pa2,
                     trackline=trackline,
                     wt_pasite2=wt2,
                     focus=focus,
                     domain=domain,
                     sitedepth=bam_site_cov,
                     logtrans=log,
                     prob=prob,
                     model=model,
                     addexpress=ade,
                     id_keep=id_keep,
                     highlight_sj=hl,
                     include_sj=inc,
                     intron_scale=intron_scale,
                     exon_scale=exon_scale,
                     bam_order=order_dic,
                     head_track=track,
                     sameymax=sameymax
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
