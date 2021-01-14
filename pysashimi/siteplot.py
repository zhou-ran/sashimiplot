#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/3/21 6:23 PM
__author__ = 'Zhou Ran'

import sys
import scipy.stats as sts

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pylab
from mpl_axes_aligner import align

from .Constant import COLOR
from .mRNA import mRNA
from .sitedepth import SiteDepth


def plot_density_single_site(read_depth_object,
                             samplename,
                             graphcoords,
                             graphToGene,
                             avx,
                             xlabel,
                             color='r',
                             nxticks=4,
                             font_size=8,
                             logtrans=None
                             ):
    """
    :param read_depth_object:
    :param mRNAs:
    :param samplename:
    :param graphcoords:
    :param graphToGene:
    :param avx:
    :param xlabel:
    :param sjthread:
    :param color:
    :param ymax:
    :param number_junctions:
    :param resolution:
    :param nxticks:
    :param font_size:
    :param numbering_font_size:
    :param junction_log_base:
    :return:
    """
    plus = read_depth_object.plus
    minus = read_depth_object.minus

    if logtrans == 'log2':
        plus = np.log2(plus + 1)
        minus = -np.log2(abs(minus) + 1)
    elif logtrans == 'log10':
        plus = np.log10(plus + 1)
        minus = -np.log10(abs(minus) + 1)
    else:
        pass

    maxheight = max(plus)
    minheight = min(minus)
    max_val = max(maxheight, abs(minheight))

    ymax = 1.1 * max_val
    ymin = 1.1 * -max_val

    pylab.fill_between(graphcoords,
                       plus,
                       y2=0,
                       color=color,
                       linewidth=0,
                       step='pre',
                       interpolate=False,
                       rasterized=True)
    pylab.fill_between(graphcoords,
                       minus,
                       y2=0,
                       color=color,
                       linewidth=0,
                       step='pre',
                       interpolate=False,
                       rasterized=True)
    for label, array_plot in zip(['plus', 'minus'], [plus, minus]):

        array_hist = np.repeat(graphcoords, np.abs(array_plot).astype(np.int))
        try:
            kde = sts.gaussian_kde(array_hist)
            fit_value = kde.pdf(graphcoords)
        except Exception:
            continue

        fit_value = fit_value / fit_value.max()
        if label == 'plus':
            plt.plot(graphcoords, fit_value * array_plot.max(), c=color, lw=1)
        else:
            plt.plot(graphcoords, fit_value * array_plot.min(), c=color, lw=1)

    # set the y limit

    avx.set_ybound(lower=ymin, upper=ymax)
    universal_yticks = pylab.linspace(ymin, ymax, 3)

    curr_yticklabels = []
    for label in universal_yticks:
        curr_yticklabels.append("{}".format(int(label)))

    avx.set_yticklabels(curr_yticklabels,
                        fontsize=font_size)
    avx.spines["left"].set_bounds(ymin, ymax)
    avx.set_yticks(universal_yticks)
    avx.yaxis.set_ticks_position('left')
    avx.spines["right"].set_visible(False)

    # ylab
    y_horz_alignment = 'right'
    avx.set_ylabel(samplename,
                   fontsize=font_size * 1.25,
                   va="center",
                   rotation="horizontal",
                   ha=y_horz_alignment,
                   labelpad=10
                   )

    avx.spines['right'].set_color('none')
    avx.spines['top'].set_color('none')

    if xlabel:
        avx.xaxis.set_ticks_position('bottom')

        max_graphcoords = max(graphcoords) - 1
        pylab.xticks(pylab.linspace(0, max_graphcoords, nxticks),
                     [graphToGene[int(x)] for x in \
                      pylab.linspace(0, max_graphcoords, nxticks)],
                     fontsize=font_size)
    else:
        avx.spines['bottom'].set_color('none')
        pylab.xticks([])

    # Here to plot the highlight site, for example pasite.
    pylab.xlim(0, max(graphcoords))


    return avx


def getScaling(tx_start,
               tx_end,
               exon_starts,
               exon_ends,
               intron_scale,
               exon_scale
               ):
    """
    Compute the scaling factor across various genic regions.
    """
    exoncoords = pylab.zeros((tx_end - tx_start + 1))

    for i in range(len(exon_starts)):
        '''
        1.22 add a if-else to solve these condition that exons were greater than the given region
        '''
        leftsite = exon_starts[i] - tx_start if exon_starts[i] - tx_start > 0 else 0
        rightsite = exon_ends[i] - tx_start if exon_ends[i] - tx_end < 0 else tx_start - tx_end

        exoncoords[leftsite: rightsite] = 1

    graphToGene = {}
    graphcoords = pylab.zeros((tx_end - tx_start + 1), dtype='f')
    x = 0

    for i in range(tx_end - tx_start + 1):
        graphcoords[i] = x
        graphToGene[int(x)] = i + tx_start
        if exoncoords[i] == 1:
            x += 1. / exon_scale
        else:
            x += 1. / intron_scale

    return graphcoords, graphToGene


def plot_density_site(read_depth_object,
                      mRNAobject,
                      fileout=None,
                      wide=8,
                      height=12,
                      logtrans=None
                      ):
    """

    :param read_depth_object:
    :param mRNAobject:
    :param strand:
    :param fileout:
    :param pasite:
    :return:
    """

    from .plot import plot_mRNAs

    txstart = mRNAobject.tstart
    txend = mRNAobject.tend

    exon_starts = mRNAobject.exonstarts
    exon_ends = mRNAobject.exonend

    intron_scale = 30
    exon_scale = 1

    plt.rcParams["figure.figsize"] = (wide, height)

    graphcoords, graphToGene = getScaling(txstart,
                                          txend,
                                          exon_starts,
                                          exon_ends,
                                          intron_scale,
                                          exon_scale
                                          )

    nfile = len(read_depth_object)
    mRNAnum = len(mRNAobject.txlst) * 2

    gs = gridspec.GridSpec(int(nfile + mRNAnum),
                           1,
                           height_ratios=[4] * nfile + [1] * mRNAnum
                           )

    for fileindex, bam_label in enumerate(read_depth_object):
        axvar = pylab.subplot(gs[fileindex, :])

        bamread = read_depth_object[bam_label]
        bamname = bam_label
        xlabel = True if fileindex == nfile - 1 else False

        try:
            color = COLOR[fileindex]
        except IndexError:
            color = COLOR[fileindex % 11]

        plot_density_single_site(bamread,
                                 bamname,
                                 graphcoords,
                                 graphToGene,
                                 axvar,
                                 xlabel,
                                 color,
                                 nxticks=4,
                                 font_size=6,
                                 logtrans=logtrans
                                 )

    pylab.subplot(gs[nfile:, :])

    plot_mRNAs(txstart,
               txend,
               mRNAobject.txlst,
               graphcoords,
               domain=False)
    if fileout:
        plt.savefig(fileout,
                    bbox_inches='tight')


def main(file):
    read_depth_object = SiteDepth.determine_depth(
        '/Volumes/bu15191450186/zr/singlecell/Microwell/SRR6954503/SRR6954503.merge.refFlag.all.q10.bam', '1', 4776354,
        4776810, "RF")
    mRNAobject = mRNA('1', 4776354, 4776810, file)

    read_depth_object = [{"BM1": read_depth_object}]

    plot_density(read_depth_object, mRNAobject, 'sashimi.color.11.pdf')
    # plt.savefig('sashimi.color.11.pdf', bbox_inches='tight')


if __name__ == '__main__':
    main(sys.argv[1])
