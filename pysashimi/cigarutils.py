#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/2/22 11:49 AM
__author__ = 'Zhou Ran'

"""
This module was used to manipulate the CIGAR string to convert the human-read information.

All relative information were compared to the reference.

NOTE:

M	BAM_CMATCH	    0
I	BAM_CINS	    1
D	BAM_CDEL	    2
N	BAM_CREF_SKIP	3
S	BAM_CSOFT_CLIP	4
H	BAM_CHARD_CLIP	5
P	BAM_CPAD	    6
=	BAM_CEQUAL	    7
X	BAM_CDIFF	    8
B	BAM_CBACK	    9

"""


def fetch_exon(start, cigar):
    """
    return the exon information, the start site must be 0-based
    :param chrom:
    :param start:
    :param cigar:
    :return:
    """
    exonbound = []

    for c, l in cigar:

        if c == 0:  # for match

            exonbound.append((start, start + l))
            start += l

        elif c == 1:  # for insert
            continue

        elif c == 2:  # for del
            start += l

        elif c == 3:  # for intron
            start += l

        elif c == 4:  # soft clip
            start += l

        else:
            continue
    return exonbound


def fetch_intron(start, cigar):
    """
    To retrieve the 'N' cigar
    :param start:
    :param cigar:
    :return:
    """
    intronbound = []

    for c, l in cigar:

        if c == 3:
            intronbound.append((start + 1, start + l))
            start += l

        elif c == 1:
            continue

        elif c == 2:
            start += l

        elif c == 4:
            start += l

        elif c == 0:
            start += l

        else:
            continue

    return intronbound
