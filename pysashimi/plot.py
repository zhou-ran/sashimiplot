#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/1/10 4:18 PM
__author__ = 'Zhou Ran'
import sys
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pylab
from loguru import logger
from matplotlib.patches import PathPatch
from matplotlib.path import Path
from scipy import stats

from .DomainCds import calculateinterval
from .Constant import COLOR
from .Constant import DOMAINFILTER
from .siteplot import plot_density_single_site
from .EM import EM
# from .color import darken_rgb

plt.switch_backend('agg')

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams["font.family"] = 'Arial'


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

    :param pts:
    :param t:
    :return:
    """
    p0, p1, p2, p3 = pts
    p0 = pylab.array(p0)
    p1 = pylab.array(p1)
    p2 = pylab.array(p2)
    p3 = pylab.array(p3)
    return p0 * (1 - t) ** 3 + 3 * t * p1 * (1 - t) ** 2 + \
        3 * t ** 2 * (1 - t) * p2 + t ** 3 * p3


def plot_density_single(read_depth_object,
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
                        highlight=None,
                        include_sj=None,
                        pasite=None,
                        wt_pasite=None,
                        pasite2=None,
                        wt_pasite2=None,
                        logtrans=None
                        ):
    """

    :param read_depth_object:
    :param samplename:
    :param graphcoords:
    :param graphToGene:
    :param avx:
    :param xlabel:
    :param sjthread:
    :param color:
    :param ymax:
    :param number_junctions:
    :param nxticks:
    :param font_size:
    :param numbering_font_size:
    :param junction_log_base:
    :param highlight:
    :param include_sj:
    :param pasite:
    :param wt_pasite:
    :param pasite2:
    :param wt_pasite2:
    :param logtrans:
    :return:
    """

    if include_sj:
        assert isinstance(
            include_sj, list), "The including junction was not a set, pls check your including input"

    if include_sj and highlight:
        if not include_sj & highlight:
            raise Exception(
                "There was no overlap bewteen including and highlight junction!")

    if highlight:
        assert isinstance(
            highlight, list), "The highlight junction was not a set, pls check your highlight input"
    else:
        highlight = set("")

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
        """
        2020.2.13 for scale mode
        """
        if maxheight == 100:
            ymax = 100
        else:
            ymax = 1.1 * maxheight
    else:
        ymax = ymax
    ymin = -.6 * ymax

    if isinstance(color, set):
        color = list(color)[0]
    else:
        pass

    pylab.fill_between(graphcoords,
                       wiggle,
                       y2=0,
                       color=color,
                       linewidth=0,
                       step='pre',
                       interpolate=False,
                       rasterized=True)

    u'''
    1.13 sorted the list by the location
    '''
    jxns_sorted_list = sorted(
        jxns.keys(), key=lambda x: int(x.split(':')[1].split('-')[0]))
    current_height = -3 * ymin / 4

    jxns_new = {}

    for sj_id, sj_counts in jxns.items():
        sj_id = sj_id.split(":")[-1].replace("-", ":")
        jxns_new[sj_id] = sj_counts
    if include_sj:
        if len(include_sj) == 2:
            psi_numerator = jxns_new.get(include_sj[0]) if jxns_new.get(
                include_sj[0]) is not None else 0
            psi_denominator = [jxns_new.get(key) for key in include_sj]
            psi = psi_numerator / \
                np.sum([0 if i is None else i for i in psi_denominator])

        elif len(include_sj) >= 3:
            psi_numerator = [jxns_new.get(key) for key in include_sj[:2]]
            psi_denominator = [jxns_new.get(key) for key in include_sj]
            psi = np.sum([0 if i is None else i for i in psi_numerator]) / np.sum(
                [0 if i is None else i for i in psi_denominator])
        else:
            psi = None

        if psi:
            psi = np.round(psi, 3)
        else:
            psi = "NA"

    # for plotted_count, jxn in enumerate(jxns_sorted_list):
    # Ensure normal sj display when including_sj mode was on

    plotted_count = 0
    for jxn in jxns_sorted_list:
        leftss, rightss = map(int, jxn.split(":")[1].split("-"))
        if round(jxns[jxn], 2) <= sjthread:
            continue
        u'''
        1.13 add a junction in half, delete this?
        just skip the junction plot
        '''
        junc_label = "{}:{}".format(leftss,
                                    rightss)

        if include_sj and junc_label not in include_sj:
            continue

        leftstatus = False
        rightstatus = False

        try:

            ss1, ss2 = [graphcoords[leftss - tx_start - 1],
                        graphcoords[rightss - tx_start + 1]]
        except IndexError:
            continue

        # draw junction on bottom

        if plotted_count % 2 == 1:
            pts = [(ss1, 0), (ss1, -current_height),
                   (ss2, -current_height), (ss2, 0)]
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

        plotted_count += 1
        if number_junctions:
            t = pylab.text(midpt[0], midpt[1], '{0}'.format(round(jxns[jxn], 2)), fontsize=numbering_font_size,
                           ha='center',
                           va='center',
                           backgroundcolor='w')

            t.set_bbox(dict(alpha=0, ))

        a = Path(pts, [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
        # here to add the highlight information for sj
        if junc_label in set(highlight):
            p = PathPatch(a, color="#282828",
                          lw=pylab.log(jxns[jxn] + 1) /
                          pylab.log(junction_log_base) * 0.5,
                          fc='none')
        else:
            p = PathPatch(a, color=color,
                          lw=pylab.log(jxns[jxn] + 1) / pylab.log(junction_log_base) * 0.5, fc='none')

        avx.add_patch(p)

    # set the y limit
    if ymax:
        max_used_yval = ymax
    else:
        max_used_yval = avx.get_ylim()[1]

    fake_ymin = -0.5 * max_used_yval
    if ymax == 100:
        avx.set_ybound(lower=fake_ymin, upper=max_used_yval)
    else:
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
    if include_sj:
        avx.set_ylabel('{}\npsi:{}'.format(samplename, psi),
                       fontsize=font_size * 1.25,
                       va="center",
                       rotation="horizontal",
                       ha=y_horz_alignment,
                       labelpad=10
                       )
    else:
        avx.set_ylabel('{}'.format(samplename),
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
                     [graphToGene[int(x)] for x in pylab.linspace(
                         0, max_graphcoords, nxticks)],
                     fontsize=font_size)

    else:
        avx.spines['bottom'].set_color('none')
        pylab.xticks([])

    # Here to plot the highlight site, for example pasite.
    if pasite:
        pasite = list(map(int, pasite.split(',')))
        wt_pasite = list(map(str, wt_pasite.split(','))
                         ) if wt_pasite else wt_pasite

        for site_index, subsite in enumerate(pasite):
            pylab.axvline(graphcoords[subsite - tx_start],
                          color="blue",
                          lw=.5)
            if wt_pasite:
                if len(pasite) != len(wt_pasite):
                    logger.error(
                        "The length of pasite and wt_pasite was not same!")
                    sys.exit()

                pylab.text(x=graphcoords[subsite - tx_start],
                           # for distinguish the overlaping weight
                           y=max_used_yval * 2 / 3,
                           s=wt_pasite[site_index],
                           color='blue',
                           fontsize=font_size)

    if pasite2:
        pasite2 = list(map(int, pasite2.split(',')))
        wt_pasite2 = list(map(str, wt_pasite2.split(','))
                          ) if wt_pasite2 else wt_pasite2

        for site_index, subsite in enumerate(pasite2):
            pylab.axvline(graphcoords[subsite - tx_start],
                          color="red",
                          lw=.5)
            if wt_pasite2:
                if len(pasite2) != len(wt_pasite2):
                    logger.error(
                        "The length of pasite and wt_pasite was not same!")
                    sys.exit()

                pylab.text(x=graphcoords[subsite - tx_start],
                           # for distinguish the overlaping weight
                           y=max_used_yval * 1 / 3,
                           s=wt_pasite2[site_index],
                           color='red',
                           fontsize=font_size)

    pylab.xlim(0, max(graphcoords))
    return avx


# def getScaling(tx_start,
#                tx_end,
#                exon_starts,
#                exon_ends,
#                intron_scale,
#                exon_scale
#                ):
#     """
#     Compute the scaling factor across various genic regions.
#     """
#     exoncoords = pylab.zeros((tx_end - tx_start + 1))
#
#     # for i in range(len(exon_starts)):
#     #     '''
#     #     1.22 add a if-else to solve these condition that exons were greater than the given region
#     #     '''
#     #     leftsite = exon_starts[i] - tx_start if exon_starts[i] - tx_start > 0 else 0
#     #     rightsite = exon_ends[i] - tx_start if exon_ends[i] - tx_end < 0 else tx_start - tx_end
#
#     #     exoncoords[leftsite: rightsite] = 1
#
#     graphToGene = {}
#     graphcoords = pylab.zeros((tx_end - tx_start + 1), dtype='f')
#     x = 0
#
#     for i in range(tx_end - tx_start + 1):
#         graphcoords[i] = x
#         graphToGene[int(x)] = i + tx_start
#         if exoncoords[i] == 1:
#             x += 1. / exon_scale
#         else:
#             x += 1. / intron_scale
#
#     return graphcoords, graphToGene

def getScaling(tx_start,
               tx_end,
               exon_starts,
               exon_ends,
               intron_scale=15,
               exon_scale=1):
    """

    :param tx_start:
    :param tx_end:
    :param exon_starts:
    :param exon_ends:
    :param intron_scale:
    :param exon_scale:
    :return:
    """
    exoncoords = np.zeros((tx_end - tx_start + 1))
    for i in range(len(exon_starts)):
        exoncoords[exon_starts[i] - tx_start: exon_ends[i] - tx_start] = 1

    graphToGene = {}
    graphcoords = np.zeros((tx_end - tx_start + 1), dtype='float64')
    # graphcoords = np.repeat(np.NaN, (tx_end - tx_start + 1))
    x = 0
    for i in range(tx_end - tx_start + 1):
        graphcoords[i] = x
        graphToGene[np.int(x)] = i + tx_start
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
    :param tx_end:
    :param graphcoords:
    :param exonwidth:
    :param yloc:
    :param xaxisloc:
    :param tid:
    :return:
    """

    RGB_tuples = ["#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
                  "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"]
    for index, subdomain in enumerate(region):
        logger.debug("plotting domain")

        if index % 2 == 0:
            u'''
            2020 4.18 to distinguish the different domain information
            '''
            y_axis_offset = 0.3
            yloc_offset_at_y = 0.05 + yloc
        else:
            y_axis_offset = -0.6
            yloc_offset_at_y = -0.05 + yloc

        dregion, domainname = subdomain
        logger.debug("{}".format(domainname))

        domain_detail, domainname = domainname.split(';;')

        '''
        1.23 add domain information, retrieve from the nature communication
        '''
        if domainname not in DOMAINFILTER:
            continue

        dregion = calculateinterval(dregion, (tx_start, tx_end))
        if not dregion:
            continue
        minsite = min(map(lambda x: x[0], dregion))
        maxsite = max(map(lambda x: x[1], dregion))
        min_ = graphcoords[0 if minsite - tx_start < 0 else minsite - tx_start]
        max_ = graphcoords[tx_end - tx_start if maxsite -
                           tx_end > 0 else maxsite - tx_start]

        for s, e in dregion:
            s = s - tx_start
            e = e - tx_start
            x = [graphcoords[s], graphcoords[e],
                 graphcoords[e], graphcoords[s]]

            y = [yloc_offset_at_y - exonwidth / 2, yloc_offset_at_y - exonwidth / 2,
                 yloc_offset_at_y + exonwidth / 2, yloc_offset_at_y + exonwidth / 2]

            '''
            1.16 domain numbers more than color numbers
            '''
            try:
                color_use = RGB_tuples[index]
            except IndexError:
                color_use = RGB_tuples[index % 8]

            pylab.fill(x, y, color_use, lw=.5, zorder=20)

        pylab.plot([min_, max_], [yloc_offset_at_y,
                                  yloc_offset_at_y], color=color_use, lw=0.2)
        pylab.text(graphcoords[s], yloc + y_axis_offset,
                   domain_detail,
                   color=color_use,
                   fontdict={'fontsize': 4},
                   bbox=dict(facecolor='none',
                             edgecolor='none')
                   )

    plt.text(xaxisloc, yloc - 0.05, '{}_domain'.format(tid), fontsize=8)


def plot_mRNAs(tx_start,
               tx_end,
               mRNAs,
               graphcoords,
               domain=True,
               focus=None,
               trackline=None,
               weight_color=None,
               id_keep=None,
               head_track=None
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
            if id_keep is not None:
                if tid not in id_keep:
                    continue
                else:
                    pass
            else:
                pass

            strand = info['strand']

            logger.debug('ploting {}'.format(tid))

            '''
            1.21 add plot the splicing plot
            '''

            toplot = ['domain', 'CDS', 'exon'] if domain else ['exon']

            for type_ in toplot:
                # minsite = min(map(lambda x: x[0], info[type_]))
                # maxsite = max(map(lambda x: x[1], info[type_]))

                logger.debug("plotting the {}".format(type_))
                try:
                    region = info[type_]
                except KeyError:
                    logger.warning(f'There was no {type_} information')
                    continue

                if not region:
                    continue
                if type_ == 'exon':
                    minsite, maxsite = info['maxinfo'][0]
                else:
                    minsite = min(map(lambda x: x[0], region))
                    maxsite = max(map(lambda x: x[1], region))

                if type_ == 'domain':
                    plotdomain(region,
                               tx_start,
                               tx_end,
                               graphcoords,
                               exonwidth,
                               yloc,
                               xaxisloc,
                               tid)
                    yloc += 1
                    continue

                region = calculateinterval(region, (tx_start, tx_end))

                if not region:
                    continue

                min_ = graphcoords[0 if minsite -
                                   tx_start < 0 else minsite - tx_start]
                max_ = graphcoords[tx_end - tx_start if maxsite -
                                   tx_end > 0 else maxsite - tx_start]

                for s, e in region:
                    s = s - tx_start
                    e = e - tx_start
                    x = [graphcoords[s], graphcoords[e],
                         graphcoords[e], graphcoords[s]]
                    y = [yloc - exonwidth / 2, yloc - exonwidth / 2,
                         yloc + exonwidth / 2, yloc + exonwidth / 2]

                    pylab.fill(x, y, '#000000', lw=.5, zorder=20)
                    pylab.plot([min_, max_], [yloc, yloc],
                               color='#000000', lw=0.5)

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
                               '|'.join(
                                   [info['symbol'], '_'.join([tid, type_])]),
                               fontsize=8)
                else:
                    xaxisloc = -1 * max(graphcoords) * 0.3
                    symbol_id = info['symbol'] if info['symbol'] else 'NAs'
                    pylab.text(xaxisloc, yloc - 0.05,
                               '|'.join([symbol_id, tid, strand]),
                               fontsize=8)
                yloc += 1
                logger.debug('{} done'.format(tid))
    if head_track:
        for label, region in head_track.items():
            if not region:
                continue
            else:
                region = calculateinterval(region, (tx_start, tx_end))

                if not region:
                    continue

                minsite = min(map(lambda x: x[0], region))
                maxsite = max(map(lambda x: x[1], region))

                min_ = graphcoords[0 if minsite -
                                    tx_start < 0 else minsite - tx_start]
                max_ = graphcoords[tx_end - tx_start if maxsite -
                                    tx_end > 0 else maxsite - tx_start]

                for s, e in region:
                    s = s - tx_start
                    e = e - tx_start
                    x = [
                        graphcoords[s],
                        graphcoords[e],
                        graphcoords[e],
                        graphcoords[s]
                        ]
                    y = [
                        yloc - exonwidth / 6,
                        yloc - exonwidth / 6,
                        yloc + exonwidth / 6,
                        yloc + exonwidth / 6
                        ]
                    pylab.fill(x, y, '#000000', lw=.5, zorder=20)

            xaxisloc = -1 * max(graphcoords) * 0.3
            pylab.text(
                xaxisloc, yloc - 0.05,
                label,
                fontsize=8
                )
            yloc += 1

    pylab.xlim(0, max(graphcoords))
    pylab.ylim(-.5, yloc + .5)
    pylab.box(on=False)

    cmap = plt.get_cmap('jet')
    color_ramp = cmap(np.linspace(0.1, 1, 10))

    if focus:
        if weight_color:
            weight_color = list(map(float, weight_color.split(',')))

        for index_, focus_ in enumerate(focus):
            l, r = list(map(int, focus_))
            fill_x = [graphcoords[l - tx_start], graphcoords[r - tx_start],
                      graphcoords[r - tx_start], graphcoords[l - tx_start]]
            fill_y = [0, 0, yloc, yloc]
            if weight_color:
                color_for_weight = color_ramp[int(10 * weight_color[index_])]
            else:
                color_for_weight = 'grey'
            pylab.fill(fill_x, fill_y, alpha=0.1, color=color_for_weight)

    if trackline:
        trackline = list(map(int, trackline.split(',')))

        for site_index, subsite in enumerate(trackline):
            pylab.axvline(graphcoords[subsite - tx_start],
                          color="blue",
                          lw=.5)
    pylab.xticks([])
    pylab.yticks([])


def plot_density(read_depth_object,
                 mRNAobject,
                 fileout,
                 sjthread=1,
                 width=8.0,
                 height=12.0,
                 colors=None,
                 pasite=None,
                 wt_pasite=None,
                 pasite2=None,
                 trackline=None,
                 wt_pasite2=None,
                 focus=None,
                 domain=True,
                 sitedepth=None,
                 logtrans=None,
                 prob=None,
                 id_keep=None,
                 bam_order=None,
                 head_track=None,
                 model=None,
                 addexpress=None,
                 sameymax=False,
                 highlight_sj=None,
                 include_sj=None,
                 intron_scale=15,
                 exon_scale=1
                 ):
    """

    :param read_depth_object:
    :param mRNAobject:
    :param fileout:
    :param sjthread:
    :param width:
    :param height:
    :param pasite:
    :param focus:
    :param domain:
    :param sitedepth:
    :param logtrans:
    :param prob:
    :param id_keep:
    :param model:
    :param addexpress:
    :param hL:
    :return:
    """
    txstart = mRNAobject.tstart
    txend = mRNAobject.tend

    # mRNAlist = mRNAobject.exon
    exon_starts = mRNAobject.exonstarts
    exon_ends = mRNAobject.exonend

    plt.rcParams["figure.figsize"] = (width, height)

    # print(txstart, txend)
    graphcoords, graphToGene = getScaling(txstart,
                                          txend,
                                          exon_starts,
                                          exon_ends,
                                          intron_scale=intron_scale,
                                          exon_scale=exon_scale
                                          )
    if sitedepth:
        nfile = 2 * len(read_depth_object)
    else:
        nfile = len(read_depth_object)

    if prob:
        nfile += 2

    # if add id_keep parameters
    if id_keep is not None:
        mRNAnum = len(id_keep) * 3
    else:
        mRNAnum = len(mRNAobject.txlst) * \
            2 if len(mRNAobject.txlst) != 0 else 1

    gs = gridspec.GridSpec(int(nfile + mRNAnum),
                           1,
                           height_ratios=[4] * nfile + [1] * mRNAnum
                           )
    if bam_order:
        colors = bam_order
        bam_labels = list(bam_order.keys())
        whether_exit = list(
            map(lambda x: x in read_depth_object.keys(), bam_labels))
        if not all(whether_exit):
            return_val = []
            for index, bool_val in enumerate(whether_exit):
                if not bool_val:
                    return_val.append(bam_labels[index])

            logger.error(
                f'{return_val}, bam order files were not all in read_depth_object, please check.')
            sys.exit(1)

    else:
        bam_labels = list(read_depth_object.keys())

    # here to scale all track into one ymax value
    if sameymax:
        ymax_val = []
        for _, bam_obj in read_depth_object.items():
            ymax_val.append(np.max(bam_obj.wiggle))
        ymax_val = np.round(max(ymax_val)) + 1
    else:
        ymax_val = None

    for file_index, bam_label in enumerate(bam_labels):
        if not sitedepth:
            x_label = True if file_index == nfile - 1 else False
            fileindex_grid = file_index
        else:
            fileindex_grid = 2 * file_index
            x_label = False

        axvar = pylab.subplot(gs[fileindex_grid, :])

        bamread = read_depth_object[bam_label]
        bamname = bam_label

        if colors:
            color = colors[bamname]
        else:
            try:
                color = COLOR[file_index]
            except IndexError:
                color = COLOR[file_index % 11]

        plot_density_single(bamread,
                            # mRNAlist,
                            bamname,
                            graphcoords,
                            graphToGene,
                            axvar,
                            x_label,
                            sjthread,
                            color,
                            ymax=ymax_val,
                            number_junctions=True,
                            nxticks=4,
                            font_size=6,
                            numbering_font_size=6,
                            junction_log_base=10,
                            pasite=pasite,
                            wt_pasite=wt_pasite,
                            pasite2=pasite2,
                            wt_pasite2=wt_pasite2,
                            logtrans=logtrans,
                            highlight=highlight_sj,
                            include_sj=include_sj
                            )
        if sitedepth:
            axvar = pylab.subplot(gs[fileindex_grid + 1, :])
            bamread = sitedepth[bam_label]
            bamname = bam_label
            x_label = True if file_index == nfile / 2 - 1 else False
            plot_density_single_site(bamread,
                                     bamname,
                                     graphcoords,
                                     graphToGene,
                                     axvar,
                                     x_label,
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
            expInfo[int(peaks) - txstart: int(peake) -
                    txstart + 1] = emMAP.pi[::-1]

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

    if wt_pasite is not None or wt_pasite2 is not None:
        weight_color = wt_pasite if wt_pasite is not None else wt_pasite2
    else:
        weight_color = None
    plot_mRNAs(txstart,
               txend,
               mRNAobject.txlst,
               graphcoords,
               domain,
               focus=focus,
               trackline=trackline,
               weight_color=weight_color,
               id_keep=id_keep,
               head_track=head_track)

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
                # figsize = c(float(width),
                #             float(height)),
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
