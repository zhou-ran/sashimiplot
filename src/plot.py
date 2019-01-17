#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/1/10 4:18 PM
__author__ = 'Zhou Ran'
import sys
import pylab
import matplotlib.gridspec as gridspec
from matplotlib.patches import PathPatch
from matplotlib.path import Path
import math
import matplotlib.pyplot as plt

plt.switch_backend('agg')
from .DomainCds import calculateinterval
from .Constant import COLOR

import logging

from logging import handlers

logger = logging.getLogger()
handler = logging.StreamHandler()
fh = handlers.RotatingFileHandler(
    'plot.log',
    mode='a+',
    encoding="utf-8",
    maxBytes=5 * 1024 * 1024,
    backupCount=2,
    delay=0
)
fh.setLevel(logging.DEBUG)
formatter = logging.Formatter(
    '%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.addHandler(fh)
logger.setLevel(logging.DEBUG)


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
                        strand,
                        samplename,
                        graphcoords,
                        graphToGene,
                        avx,
                        xlabel,
                        color='r',
                        ymax=None,
                        number_junctions=True,
                        resolution=.5,
                        nxticks=4,
                        font_size=8,
                        numbering_font_size=6,
                        junction_log_base=10):
    # extract data from read_depth_object 4702888, 4784871

    tx_start = read_depth_object.low
    tx_end = read_depth_object.high
    chrom = read_depth_object.chrm
    wiggle = read_depth_object.wiggle
    jxns = read_depth_object.junctions_dict

    maxheight = max(wiggle)

    if ymax is None:
        ymax = 1.1 * maxheight
    else:
        ymax = ymax
    ymin = -.6 * ymax

    # Reduce memory footprint by using incremented graphcoords.
    compressed_x = []
    compressed_wiggle = []
    prevx = graphcoords[0]
    tmpval = []

    for i in range(len(graphcoords)):
        tmpval.append(wiggle[i])
        if abs(graphcoords[i] - prevx) > resolution:
            compressed_wiggle.append(pylab.mean(tmpval))
            compressed_x.append(prevx)
            prevx = graphcoords[i]
            tmpval = []

    pylab.fill_between(compressed_x, compressed_wiggle, \
                       y2=0, color=color, lw=0)

    sslists = []
    for mRNA in mRNAs:
        tmp = []
        for s, e in mRNA:
            tmp.extend([s, e])
        sslists.append(tmp)

    # sort the junctions by intron length for better plotting look
    # jxns_sorted_list = sorted(jxns.keys(),cmp=junc_comp_function)
    # ['1:4774516-4777525', '1:4774516-4776410', '1:4777648-4782568', '1:4776704-4903966', '1:4776801-4777525',
    #  '1:4782060-4782568', '1:4782301-4782568', '1:4782733-4783951', '1:4782733-4785573', '1:4782733-4784078',
    #  '1:4784105-4785573', '1:4784105-4785618', '1:4784100-4785573']

    u'''
    1.13 sorted the list by the location
    '''
    jxns_sorted_list = sorted(jxns.keys(), key=lambda x: int(x.split(':')[1].split('-')[0]))
    current_height = -3 * ymin / 4

    for plotted_count, jxn in enumerate(jxns_sorted_list):
        leftss, rightss = map(int, jxn.split(":")[1].split("-"))

        u'''
        1.13 add a junction in half
        '''
        leftstatus = False
        rightstatus = False
        try:
            ss1, ss2 = [graphcoords[leftss - tx_start - 1], graphcoords[rightss - tx_start]]
        except IndexError:

            leftsite = leftss - tx_start
            rightsite = rightss - tx_end

            if leftsite < 0 and rightsite > 0:
                continue
            else:
                if leftsite < 0:
                    ss1, ss2 = [graphcoords[0], graphcoords[rightss - tx_start]]
                    leftstatus = True
                else:
                    ss1, ss2 = [graphcoords[leftss - tx_start - 1], graphcoords[tx_end - tx_start]]
                    rightstatus = True

        mid = (ss1 + ss2) / 2

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
                rightdens = wiggle[rightss - tx_start]
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
                           va='center', backgroundcolor='w')
            t.set_bbox(dict(alpha=0, ))

        a = Path(pts, [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
        p = PathPatch(a, ec=color, lw=pylab.log(jxns[jxn] + 1) / pylab.log(junction_log_base) * 0.5, fc='none')
        avx.add_patch(p)

    # set the y limit
    max_used_yval = avx.get_ylim()[1]
    fake_ymin = -0.5 * max_used_yval

    avx.set_ybound(lower=fake_ymin, upper=1.2 * max_used_yval)
    universal_yticks = pylab.linspace(0, max_used_yval, 2 + 1)
    universal_ticks = map(math.ceil, universal_yticks)
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
            # if label % 1 != 0:
            #     curr_yticklabels.append("%.1f" % (label))
            # else:
            #     curr_yticklabels.append("%d" % (label))
    avx.set_yticklabels(curr_yticklabels,
                        fontsize=font_size)
    avx.spines["left"].set_bounds(0, max_used_yval)
    avx.set_yticks(universal_yticks)
    avx.yaxis.set_ticks_position('left')
    avx.spines["right"].set_color('none')

    # ylab
    y_horz_alignment = 'right'
    avx.set_ylabel(samplename,
                   fontsize=font_size,
                   va="center",
                   rotation="horizontal",
                   ha=y_horz_alignment,
                   labelpad=10)

    # Format plot
    # ylim(ymin, ymax)
    # axvar.spines['left'].set_bounds(0, ymax)

    avx.spines['right'].set_color('none')
    avx.spines['top'].set_color('none')

    if xlabel:
        avx.xaxis.set_ticks_position('bottom')
        pylab.xlabel('Genomic coordinate (%s), "%s" strand' % (chrom,
                                                               strand),
                     fontsize=font_size)
        max_graphcoords = max(graphcoords) - 1
        pylab.xticks(pylab.linspace(0, max_graphcoords, nxticks),
                     [graphToGene[int(x)] for x in \
                      pylab.linspace(0, max_graphcoords, nxticks)],
                     fontsize=font_size)
    else:
        avx.spines['bottom'].set_color('none')
        pylab.xticks([])

    pylab.xlim(0, max(graphcoords))
    return avx


def getScaling(tx_start,
               tx_end,
               strand,
               exon_starts,
               exon_ends,
               intron_scale,
               exon_scale,
               reverse_minus):
    """
    Compute the scaling factor across various genic regions.
    """
    exoncoords = pylab.zeros((tx_end - tx_start + 1))
    for i in range(len(exon_starts)):
        exoncoords[exon_starts[i] - tx_start: exon_ends[i] - tx_start] = 1

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

    # if strand == '+' or not reverse_minus:
    #     for i in range(tx_end - tx_start + 1):
    #         graphcoords[i] = x
    #         graphToGene[int(x)] = i + tx_start
    #         if exoncoords[i] == 1:
    #             x += 1. / exon_scale
    #         else:
    #             x += 1. / intron_scale
    # else:
    #     for i in range(tx_end - tx_start + 1):
    #         graphcoords[-(i + 1)] = x
    #         graphToGene[int(x)] = tx_end - i + 1
    #         if exoncoords[-(i + 1)] == 1:
    #             x += 1. / exon_scale
    #         else:
    #             x += 1. / intron_scale
    return graphcoords, graphToGene


def plotdomain(region,
               tx_start,
               tx_end,
               graphcoords,
               exonwidth,
               yloc,
               xaxisloc,
               tid,
               strand):
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
        logging.debug("ploting domain")
        dregion, domainname = subdomain
        logging.debug("{}".format(domainname))
        print(dregion)
        if "chain" in domainname: continue
        if "conflict" in domainname: continue
        if "variant" in domainname: continue
        try:
            dregion = calculateinterval(dregion, (tx_start, tx_end))
        except TypeError:
            print(dregion, (tx_start, tx_end))

        minsite = min(map(lambda x: x[0], dregion))
        maxsite = max(map(lambda x: x[1], dregion))
        min_ = graphcoords[0 if minsite - tx_start < 0 else minsite - tx_start]
        max_ = graphcoords[tx_end - tx_start if maxsite - tx_end > 0 else maxsite - tx_start]
        # if strand == "+":
        #     min_ = graphcoords[0 if minsite - tx_start < 0 else minsite - tx_start]
        #     max_ = graphcoords[tx_end - tx_start if maxsite - tx_end > 0 else maxsite - tx_start]
        # else:
        #     min_ = graphcoords[-1 if maxsite - tx_end > 0 else maxsite - tx_end - 1]
        #     max_ = graphcoords[0 if minsite - tx_start < 0 else minsite - tx_end - 1]

        for s, e in dregion:
            s = s - tx_start
            e = e - tx_start
            x = [graphcoords[s], graphcoords[e], graphcoords[e], graphcoords[s]]

            y = [yloc - exonwidth / 2, yloc - exonwidth / 2, \
                 yloc + exonwidth / 2, yloc + exonwidth / 2]
            # pylab.fill(x, y, RGB_tuples[index], lw=.5, zorder=20)
            # pylab.fill(x, y, 'k', lw=.5, zorder=20)
            # plot([min_ * 1.05, max_], [yloc, yloc], color='k', lw=0.5)
            # pylab.plot([min_, max_], [yloc, yloc], color='k', lw=0.5)

            '''
            1.16 domain numbers more than color numbers
            '''
            try:
                pylab.fill(x, y, RGB_tuples[index], lw=.5, zorder=20)
                pylab.plot([min_, max_], [yloc, yloc], color=RGB_tuples[index], lw=0.2)
            except IndexError:
                pylab.fill(x, y, RGB_tuples[index % 8], lw=.5, zorder=20)
                pylab.plot([min_, max_], [yloc, yloc], color=RGB_tuples[index % 8], lw=0.2)
    plt.text(xaxisloc, yloc - 0.05, '{};domain'.format(tid), fontsize=8)


def plot_mRNAs(tx_start, tx_end, mRNAs, graphcoords, reverse_minus):
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
            # print(tid)
            logging.debug('ploting {}'.format(tid))
            toplot = ['domain', 'CDS', 'exon']
            print(info)
            for type_ in toplot:
                logging.debug("ploting the {}".format(type_))
                try:
                    region = info[type_]
                except ValueError:
                    continue

                if not region:
                    continue

                if type_ == 'domain':
                    plotdomain(region, tx_start, tx_end, graphcoords, exonwidth, yloc, xaxisloc, tid, strand)
                    yloc += 1
                    continue

                region = calculateinterval(region, (tx_start, tx_end))

                minsite = min(map(lambda x: x[0], region))
                maxsite = max(map(lambda x: x[1], region))
                # if strand == "+":
                #     min_ = graphcoords[0 if minsite - tx_start < 0 else minsite - tx_start]
                #     max_ = graphcoords[tx_end - tx_start if maxsite - tx_end > 0 else maxsite - tx_start]
                # else:
                #     max_ = graphcoords[-1 if maxsite - tx_end > 0 else maxsite - tx_end - 1]
                #     min_ = graphcoords[0 if minsite - tx_start < 0 else minsite - tx_end - 1]
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

                logging.debug('{} ploting intron'.format(tid))

                spread = .5 * (max_ - min_) / narrows
                for i in range(narrows):
                    loc = float(i) * max(graphcoords) / narrows + min_
                    if loc >= max_:
                        # print(region, loc, x, min_, max_, tx_start, tx_end, graphcoords)
                        break

                    if strand == '+' or reverse_minus:
                        x = [loc - spread, loc, loc - spread]
                    else:
                        x = [loc + spread, loc, loc + spread]
                    y = [yloc - exonwidth / 5, yloc, yloc + exonwidth / 5]
                    pylab.plot(x, y, lw=.5, color='#000000')
                logging.debug('{} ploting tid'.format(tid))
                pylab.text(xaxisloc, yloc - 0.05, ';'.join([tid, type_]), fontsize=8)
                yloc += 1
                logging.debug('{} done'.format(tid))
    pylab.xlim(0, max(graphcoords))
    pylab.ylim(-.5, yloc + .5)
    pylab.box(on=False)
    pylab.xticks([])
    pylab.yticks([])


def plot_density(read_depth_object,
                 mRNAobject,
                 strand,
                 fileout,
                 pasite=None
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
    reverse_minus = False

    graphcoords, graphToGene = getScaling(txstart, txend,
                                          strand, exon_starts,
                                          exon_ends, intron_scale,
                                          exon_scale, reverse_minus)

    nfile = len(read_depth_object)
    gs = gridspec.GridSpec(int(nfile + nfile * 0.8 + 1), 1)

    for fileindex, bamfileinfo in enumerate(read_depth_object):
        axvar = pylab.subplot(gs[fileindex, :])
        bamread = list(bamfileinfo.values())[0]
        bamname = list(bamfileinfo.keys())[0]
        xlabel = True if fileindex == nfile - 1 else False

        try:
            color = COLOR[fileindex]
        except IndexError:
            color = COLOR[fileindex % 11]

        plot_density_single(bamread,
                            mRNAlist,
                            strand,
                            bamname,
                            graphcoords,
                            graphToGene,
                            axvar,
                            xlabel,
                            color,
                            ymax=None,
                            number_junctions=True,
                            resolution=.5,
                            nxticks=4,
                            font_size=6,
                            numbering_font_size=6,
                            junction_log_base=10
                            )

    pylab.subplot(gs[nfile:, :])

    plot_mRNAs(txstart,
               txend,
               mRNAobject.txlst,
               graphcoords,
               reverse_minus)

    if pasite:
        plt.axvline(graphcoords[pasite - txstart], lw=.5)
    plt.savefig(fileout, bbox_inches='tight')


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
