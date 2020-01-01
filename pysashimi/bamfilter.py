#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/3/27 9:07 PM
__author__ = 'Zhou Ran'


def pafilter(reads, peak):
    """
    polyA filter for bam file
    :param reads:
    :param peak:
    :return:
    """
    chrom, st, en, strand = peak.split(':')
    cigarinfo = reads.cigar
    sequence = reads.seq

    if reads.is_read1:
        u'''
        To check the R1 whether there were soft clip or there were polyA or polyT on sequence?
        '''

        if not reads.is_reverse and strand == '+':
            reads_site = reads.reference_start + 1

        elif reads.is_reverse and strand == '-':
            reads_site = reads.reference_end + 1

        else:
            return False

    else:
        if reads.is_reverse and strand == '-':

            reads_site = reads.reference_end + 1
            if cigarinfo[0][0] != 4:
                return False
            else:
                if 'TTTTTT' not in sequence:
                    return False
                elif cigarinfo[-1][0] == 4:
                    return False
                else:
                    pass
        elif not reads.is_reverse and strand == '+':
            reads_site = reads.reference_start + 1
            if cigarinfo[-1][0] != 4:
                return False
            else:
                if 'AAAAAAAAAAAA' not in sequence:
                    return False
                elif cigarinfo[0][0] == 4:
                    return False
                else:
                    pass
                    # return (read.reference_end + 1, '+')
                    # return (read.reference_end, '+')
        else:
            return False

    if int(st) <= reads_site <= int(en):
        return reads_site

    return False
