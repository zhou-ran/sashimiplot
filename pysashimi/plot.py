#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/1/10 4:18 PM
__author__ = 'Zhou Ran'
import sys

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pylab
from loguru import logger
from matplotlib.patches import PathPatch
from matplotlib.path import Path
from scipy import stats

plt.switch_backend('agg')

from .DomainCds import calculateinterval
from .Constant import COLOR
from .Constant import DOMAINFILTER
from .siteplot import plot_density_single_site
from .EM import EM


def r2(x, y):
    """
    calculate the r-square
    :param x:
    :param y:
    :return:
    """
    return stats.pearsonr(x, y)[0] ** 2


def cubic_bezier(pts, t):
    """
    Get points in a cubic bezier.
    """
    p0, p1, p2, p3 = pts
    p0 = pylab.array(p0)
    p1 = pylab.array(p1)
    p2 = pylab.array(p2)
    p3 = pylab.array(p3)
    return p0 * (1 - t) ** 3 + 3 * t * p1 * (1 - t) ** 2 + \
           3 * t ** 2 * (1 - t) * p2 + t ** 3 * p3


def plot_density_single(read_depth_object,
                        mRNAs,
                        samplename,
                        graphcoords,
                        graphToGene,
                        avx,
                        xlabel,
                        sjthread=1,
                        color='r',
                        ymax=None,
                        number_junctions=True,
                        nxticks=4,
                        font_size=8,
                        numbering_font_size=6,
                        junction_log_base=10,
                        pasite=None,
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

    # TODO, if given a known junction list, highligh the given junction line
    tx_start = read_depth_object.low

    wiggle = read_depth_object.wiggle
    jxns = read_depth_object.junctions_dict

    if logtrans == 'log2':
        wiggle = np.log2(wiggle + 1)
    elif logtrans == 'log10':
        wiggle = np.log10(wiggle + 1)
    else:
        pass

    maxheight = max(wiggle)

    if ymax is None:
        ymax = 1.1 * maxheight
    else:
        ymax = ymax
    ymin = -.6 * ymax

    pylab.fill_between(graphcoords,
                       wiggle,
                       y2=0,
                       color=color,
                       linewidth=0,
                       step='pre',
                       interpolate=False)

    u'''
    1.13 sorted the list by the location
    '''
    jxns_sorted_list = sorted(jxns.keys(), key=lambda x: int(x.split(':')[1].split('-')[0]))
    current_height = -3 * ymin / 4

    for plotted_count, jxn in enumerate(jxns_sorted_list):

        leftss, rightss = map(int, jxn.split(":")[1].split("-"))
        if round(jxns[jxn], 2) <= sjthread:
            continue
        u'''
        1.13 add a junction in half, delete this?
        '''

        leftstatus = False
        rightstatus = False

        try:
            ss1, ss2 = [graphcoords[leftss - tx_start - 1], graphcoords[rightss - tx_start + 1]]
        except IndexError:
            continue

        # draw junction on bottom

        if plotted_count % 2 == 1:
            pts = [(ss1, 0), (ss1, -current_height), (ss2, -current_height), (ss2, 0)]
            midpt = cubic_bezier(pts, .5)

        # draw junction on top
        else:
            if not leftstatus:
                leftdens = wiggle[leftss - tx_start - 1]
            else:
                leftdens = maxheight

            if not rightstatus:
                rightdens = wiggle[rightss - tx_start + 1]
            else:
                rightdens = maxheight

            pts = [(ss1, leftdens),
                   (ss1, leftdens + current_height),
                   (ss2, rightdens + current_height),
                   (ss2, rightdens)]
            midpt = cubic_bezier(pts, .5)

        if number_junctions:
            t = pylab.text(midpt[0], midpt[1], '{0}'.format(round(jxns[jxn], 2)), fontsize=numbering_font_size,
                           ha='center',
                           va='center',
                           backgroundcolor='w')

            t.set_bbox(dict(alpha=0, ))

        a = Path(pts, [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
        p = PathPatch(a, ec=color, lw=pylab.log(jxns[jxn] + 1) / pylab.log(junction_log_base) * 0.5, fc='none')
        avx.add_patch(p)

    # set the y limit
    max_used_yval = avx.get_ylim()[1]
    fake_ymin = -0.5 * max_used_yval

    avx.set_ybound(lower=fake_ymin, upper=1.2 * max_used_yval)
    universal_yticks = pylab.linspace(0, max_used_yval, 2 + 1)
    # universal_ticks = map(math.ceil, universal_yticks)
    curr_yticklabels = []
    for label in universal_yticks:
        if label <= 0:
            # Exclude label for 0
            curr_yticklabels.append("")
        else:
            u'''
            1.9, force int y axis 
            '''
            curr_yticklabels.append("{}".format(int(label)))

    avx.set_yticklabels(curr_yticklabels,
                        fontsize=font_size)
    avx.spines["left"].set_bounds(0, max_used_yval)
    avx.set_yticks(universal_yticks)
    avx.yaxis.set_ticks_position('left')
    avx.spines["right"].set_color('none')

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

        max_graphcoords = max(graphcoords)

        pylab.xticks(pylab.linspace(0, max_graphcoords, nxticks),
                     [graphToGene[int(x)] for x in \
                      pylab.linspace(0, max_graphcoords, nxticks)],
                     fontsize=font_size)

    else:
        avx.spines['bottom'].set_color('none')
        pylab.xticks([])

    # Here to plot the highlight site, for example pasite.
    if pasite:
        pasite = map(int, pasite.split(','))
        for subsite in pasite:
            pylab.axvline(graphcoords[subsite - tx_start], lw=.5)

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


def plotdomain(region,
               tx_start,
               tx_end,
               graphcoords,
               exonwidth,
               yloc,
               xaxisloc,
               tid
               ):
    """

    :param region:
    :param tx_start:
    :param graphcoords:
    :param exonwidth:
    :param yloc:
    :param min_:
    :param max_:
    :return:
    """

    RGB_tuples = ["#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"]
    for index, subdomain in enumerate(region):
        logger.debug("ploting domain")

        dregion, domainname = subdomain
        logger.debug("{}".format(domainname))

        domainname = domainname.split(';;')[-1]

        '''
        1.23 add domain information, retrieve from the nature communication
        '''
        if domainname not in DOMAINFILTER:
            continue

        dregion = calculateinterval(dregion, (tx_start, tx_end))

        minsite = min(map(lambda x: x[0], dregion))
        maxsite = max(map(lambda x: x[1], dregion))
        min_ = graphcoords[0 if minsite - tx_start < 0 else minsite - tx_start]
        max_ = graphcoords[tx_end - tx_start if maxsite - tx_end > 0 else maxsite - tx_start]

        for s, e in dregion:
            s = s - tx_start
            e = e - tx_start
            x = [graphcoords[s], graphcoords[e], graphcoords[e], graphcoords[s]]

            y = [yloc - exonwidth / 2, yloc - exonwidth / 2, \
                 yloc + exonwidth / 2, yloc + exonwidth / 2]

            '''
            1.16 domain numbers more than color numbers
            '''
            try:
                pylab.fill(x, y, RGB_tuples[index], lw=.5, zorder=20)
                pylab.plot([min_, max_], [yloc, yloc], color=RGB_tuples[index], lw=0.2)
            except IndexError:
                pylab.fill(x, y, RGB_tuples[index % 8], lw=.5, zorder=20)
                pylab.plot([min_, max_], [yloc, yloc], color=RGB_tuples[index % 8], lw=0.2)
    plt.text(xaxisloc, yloc - 0.05, '{}_domain'.format(tid), fontsize=8)


def plot_mRNAs(tx_start,
               tx_end,
               mRNAs,
               graphcoords,
               domain=True,
               focus=None
               ):
    """
    Draw the gene structure.
    """

    xaxisloc = -1 * max(graphcoords) * 0.4
    yloc = 0
    exonwidth = .3
    narrows = 50
    # test the size of graph
    # plt.figure(figsize=(8, 4))

    for allinfo in mRNAs:
        for tid, info in allinfo.items():
            strand = info['strand']
            # print(info)
            minsite, maxsite = info['maxinfo'][0]
            logger.debug('ploting {}'.format(tid))
            '''
            1.21 add plot the splicing plot
            '''
            toplot = ['domain', 'CDS', 'exon'] if domain else ['exon']

            for type_ in toplot:

                logger.debug("ploting the {}".format(type_))
                try:
                    region = info[type_]
                except KeyError:
                    continue

                if not region:
                    continue

                if type_ == 'domain':
                    plotdomain(region, tx_start, tx_end, graphcoords, exonwidth, yloc, xaxisloc, tid)
                    yloc += 1
                    continue

                region = calculateinterval(region, (tx_start, tx_end))

                min_ = graphcoords[0 if minsite - tx_start < 0 else minsite - tx_start]
                max_ = graphcoords[tx_end - tx_start if maxsite - tx_end > 0 else maxsite - tx_start]

                for s, e in region:
                    s = s - tx_start
                    e = e - tx_start
                    x = [graphcoords[s], graphcoords[e], graphcoords[e], graphcoords[s]]
                    y = [yloc - exonwidth / 2, yloc - exonwidth / 2, \
                         yloc + exonwidth / 2, yloc + exonwidth / 2]

                    pylab.fill(x, y, '#000000', lw=.5, zorder=20)
                    pylab.plot([min_, max_], [yloc, yloc], color='#000000', lw=0.5)

                logger.debug('{} ploting intron'.format(tid))

                spread = .5 * (max_ - min_) / narrows
                for i in range(narrows):
                    loc = float(i) * max(graphcoords) / narrows + min_
                    if loc >= max_:
                        break

                    if strand == '+':
                        x = [loc - spread, loc, loc - spread]
                    else:
                        x = [loc + spread, loc, loc + spread]
                    y = [yloc - exonwidth / 5, yloc, yloc + exonwidth / 5]
                    pylab.plot(x, y, lw=.5, color='#000000')
                logger.debug('{} ploting tid'.format(tid))

                if domain:
                    pylab.text(xaxisloc, yloc - 0.05,
                               '|'.join([info['symbol'], '_'.join([tid, type_])]),
                               fontsize=8)
                else:
                    xaxisloc = -1 * max(graphcoords) * 0.3
                    pylab.text(xaxisloc, yloc - 0.05,
                               '|'.join([info['symbol'], tid]),
                               fontsize=8)
                yloc += 1
                logger.debug('{} done'.format(tid))

    pylab.xlim(0, max(graphcoords))
    pylab.ylim(-.5, yloc + .5)
    pylab.box(on=False)
    if focus:
        for focus_ in focus:
            l, r = list(map(int, focus_))
            fill_x = [graphcoords[l - tx_start], graphcoords[r - tx_start],
                      graphcoords[r - tx_start], graphcoords[l - tx_start]]
            fill_y = [0, 0, yloc, yloc]
            pylab.fill(fill_x, fill_y, alpha=0.1, color='grey')

    pylab.xticks([])
    pylab.yticks([])


def plot_density(read_depth_object,
                 mRNAobject,
                 fileout,
                 sjthread=1,
                 wide=8,
                 height=12,
                 pasite=None,
                 focus=None,
                 domain=True,
                 sitedepth=None,
                 logtrans=None,
                 prob=None,
                 model=None,
                 addexpress=None
                 ):
    """

    :param read_depth_object:
    :param mRNAobject:
    :param strand:
    :param fileout:
    :param pasite:
    :return:
    """
    txstart = mRNAobject.tstart
    txend = mRNAobject.tend

    mRNAlist = mRNAobject.exon
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
    if sitedepth:
        nfile = 2 * len(read_depth_object)
    else:
        nfile = len(read_depth_object)

    if prob:
        nfile += 2
    mRNAnum = len(mRNAobject.txlst) * 2 if len(mRNAobject.txlst) != 0 else 1

    gs = gridspec.GridSpec(int(nfile + mRNAnum),
                           1,
                           height_ratios=[4] * nfile + [1] * mRNAnum
                           )

    for fileindex, bamfileinfo in enumerate(read_depth_object):
        if not sitedepth:
            xlabel = True if fileindex == nfile - 1 else False
            fileindex_grid = fileindex
        else:
            fileindex_grid = 2 * fileindex
            xlabel = False

        axvar = pylab.subplot(gs[fileindex_grid, :])

        bamread = list(bamfileinfo.values())[0]
        bamname = list(bamfileinfo.keys())[0]

        try:
            color = COLOR[fileindex]
        except IndexError:
            color = COLOR[fileindex % 11]

        plot_density_single(bamread,
                            mRNAlist,
                            bamname,
                            graphcoords,
                            graphToGene,
                            axvar,
                            xlabel,
                            sjthread,
                            color,
                            ymax=None,
                            number_junctions=True,
                            nxticks=4,
                            font_size=6,
                            numbering_font_size=6,
                            junction_log_base=10,
                            pasite=pasite,
                            logtrans=logtrans
                            )
        if sitedepth:
            axvar = pylab.subplot(gs[fileindex_grid + 1, :])
            bamfileinfo = sitedepth[fileindex]
            bamread = list(bamfileinfo.values())[0]
            bamname = list(bamfileinfo.keys())[0]
            xlabel = True if fileindex == nfile / 2 - 1 else False

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
    # modelrate(chrom, peak_s, peak_e, plotregion_s, plotregion_e, strand, model):
    """ For CNN
    if prob:
        chrom, peaks, peake, strand = prob.split(':')
        if addexpress == True:
            print('Adding expression')
            if strand == "+":
                genemtx = bamread.plus
            else:
                genemtx = bamread.minus
            probdat = modelrate_(
                chrom,
                int(peaks), int(peake), txstart, txend, strand,
                model, genemtx=genemtx)
        else:
            probdat = modelrate(
                chrom,
                int(peaks), int(peake), txstart, txend, strand,
                model)

        ax = pylab.subplot(gs[fileindex_grid + 2, :])

        pylab.plot(graphcoords, probdat)
        plus = bamread.plus * probdat
        minus = bamread.minus * probdat

        # plot 0 if a site was smaller than 0.5
        probdat[probdat < 0.5] = 0
        pylab.plot(graphcoords, probdat)

        ax.get_xaxis().set_ticks([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        pylab.xlim(0, max(graphcoords))
        pylab.ylim(0, 1)

        ax = pylab.subplot(gs[fileindex_grid + 3, :])
        # for weigh
        pylab.plot(graphcoords, plus)
        pylab.plot(graphcoords, minus)
        ax.get_xaxis().set_ticks([])
        # ax.get_yaxis().set_ticks([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        # ax.spines['left'].set_visible(False)
        pylab.xlim(0, max(graphcoords))
    """

    """For EM"""
    if prob:
        chrom, peaks, peake, strand = prob.split(':')
        plusCounts = bamread.plus
        minusCounts = bamread.minus
        if np.sum(plusCounts) > - np.sum(minusCounts):
            strand = '+'
            expInfo = plusCounts
        else:
            strand = '-'
            expInfo = -minusCounts[::-1]

        _expInfo = expInfo[int(peaks) - txstart: int(peake) - txstart + 1]
        emMAP = EM(_expInfo)
        emMAP.fit()
        if strand == '+':
            expInfo[int(peaks) - txstart: int(peake) - txstart + 1] = emMAP.pi
        else:
            expInfo[int(peaks) - txstart: int(peake) - txstart + 1] = emMAP.pi[::-1]

        ax = pylab.subplot(gs[fileindex_grid + 2, :])

        plt.scatter(graphcoords, expInfo, s=0.5, color='blue')
        for axvline in graphcoords[expInfo > 0.1]:
            plt.axvline(x=axvline, color='black', linewidth=.5)

        ax.get_xaxis().set_ticks([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        pylab.xlim(0, max(graphcoords))
        pylab.ylim(0, 1)

    pylab.subplot(gs[nfile:, :])

    plot_mRNAs(txstart,
               txend,
               mRNAobject.txlst,
               graphcoords,
               domain,
               focus=focus)

    # if pasite:
    #     pasite = map(int, pasite.split(','))
    #     for subsite in pasite:
    #         plt.axvline(graphcoords[subsite - txstart], lw=.5)

    '''
    1.23 add focus line into AS
    '''

    # if focus:
    #     assert len(focus) == 2, "Got a error on focus region, pls check the focus region you have given"
    #     l, r = list(map(int, focus))
    #     fill_l = [graphcoords[l - txstart], graphcoords[r - txstart],
    #               graphcoords[r - txstart], graphcoords[l - txstart]]
    #
    #     pylab.fill(fill_l)
    # plt.axvline(graphcoords[l - txstart], lw=.5)
    # plt.axvline(graphcoords[r - txstart], lw=.5)

    plt.savefig(fileout,
                bbox_inches='tight')

    """ This is the older code for analysis the cnn model
    if prob:
        fig, (ax1, ax2) = pylab.subplots(2, 1)
        sns.regplot(np.log2(abs(bamread.plus) + 1), probdat, ax=ax1)

        if sum(bamread.plus) != 0:
            rsquare = stats.pearsonr(np.log2(abs(bamread.plus) + 1), probdat)[0] ** 2
        else:
            rsquare = 'NaN'
        ax1.set_ylabel('plus(+) strand')
        ax1.set_xlabel('log2(count + 1)|R^2: {}'.format(round(rsquare, 3) if rsquare != 'NaN' else rsquare))

        sns.regplot(np.log2(abs(bamread.minus) + 1), probdat, ax=ax2)

        if sum(bamread.minus) != 0:
            rsquare = stats.pearsonr(np.log2(abs(bamread.minus) + 1), probdat)[0] ** 2
        else:
            rsquare = 'NaN'

        ax2.set_ylabel('minus(-) strand')
        ax2.set_xlabel('log2(count + 1)|R^2: {}'.format(round(rsquare, 3) if rsquare != 'NaN' else rsquare))

        fig.savefig(fileout.replace('pdf', 'cor.pdf'),
                    bbox_inches='tight')

    fig, (ax1, ax2) = pylab.subplots(2, 1)
    ax1.scatter(np.arange(0, len(bamread.plus)), np.log2(np.array(bamread.plus) + 1))
    ax2.scatter(np.arange(0, len(bamread.minus)), np.log(-np.array(bamread.minus) + 1))
    fig.savefig(fileout.replace('pdf', 'readdistribution.pdf'),
                bbox_inches='tight')
    """


def main(file):
    from ReadDepth import ReadDepth
    from mRNA import mRNA

    read_depth_object = ReadDepth.determine_depth('/Users/zhouran/opt/proj/pysashimi/test/Chr1.bam', '1', 4773206,
                                                  4785739)
    mRNAobject = mRNA('1', 4773206, 4785739, file)
    strand = '+'

    read_depth_object = [{"test1": read_depth_object},
                         {"test2": read_depth_object}]

    plot_density(read_depth_object, mRNAobject, strand, 4781892)
    plt.savefig('sashimi.color.11.pdf', bbox_inches='tight')


if __name__ == '__main__':
    main(sys.argv[1])
