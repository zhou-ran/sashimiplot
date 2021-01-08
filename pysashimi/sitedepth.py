#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/3/21 5:20 PM
__author__ = 'Zhou Ran'

"""
This class is for plotting the site class
"""

import numpy as np
import pysam

from .bamfilter import pafilter
from .utils import checkbam


class SiteDepth:
    def __init__(self,
                 chrom,
                 low,
                 high,
                 plus,
                 minus):
        self.low = low
        self.high = high
        self.chrom = chrom
        self.minus = minus
        self.plus = plus

    @classmethod
    def determine_depth(cls,
                        bam_file_path,
                        chrom,
                        start_coord,
                        end_coord,
                        libtype,
                        singlestrand=None,
                        readFilter=None):

        """
        calculate the site coverage based the libtype and the input strand.
        FR means: the R1 is reversed to the gene strand, and the R2 is same to the gene strand

        :param bam_file_path:
        :param chrom:
        :param start_coord:
        :param end_coord:
        :param libtype:
        :param strand:
        :return:
        """
        checkbam(bam_file_path)
        assert libtype in ['RF', 'FR'], "The library type must be one of RF or FR"
        if singlestrand:
            assert singlestrand in ['R1', 'R2'], "The single strand mode must be R1 or R2"

        try:
            bam_file = pysam.Samfile(bam_file_path, 'rb', ignore_truncation=True)
            chrom_id = set(map(lambda x: x['SN'], bam_file.header['SQ']))

            if chrom not in chrom_id:
                if chrom.startswith('chr'):
                    if chrom.replace('chr', '') in chrom_id:
                        chrom_id = chrom.replace('chr', '')
                    else:
                        chrom_id = chrom
                else:
                    if chrom == 'MT' and 'chrom' in chrom_id:
                        chrom_id = 'chrom'
                    elif f'chr{chrom}' in chrom_id:
                        chrom_id = f'chr{chrom}'
                    else:
                        chrom_id = chrom
            else:
                chrom_id = chrom

            relevant_reads = bam_file.fetch(reference=chrom_id, start=start_coord, end=end_coord)

            plus = np.zeros(end_coord - start_coord + 1, dtype='f')
            minus = np.zeros(end_coord - start_coord + 1, dtype='f')

            for read in relevant_reads:
                if readFilter:
                    pasite = pafilter(read, readFilter)
                    if not pasite:
                        continue

                if singlestrand == 'R1' and read.is_read2:
                    continue
                elif singlestrand == 'R2' and read.is_read1:
                    continue
                else:
                    pass

                # make sure that the read can be used
                cigar_string = read.cigar

                # each read must have a cigar string
                if cigar_string == None:
                    continue
                # if not start_coord <= read.reference_start + 1 <= end_coord and not start_coord <=
                # read.reference_end + 1 <= end_coord: continue
                if read.is_reverse:
                    if read.is_read1:
                        if not start_coord <= read.reference_end + 1 <= end_coord:
                            continue
                        plus[read.reference_end - start_coord + 1] += 1
                    else:
                        if not start_coord <= read.reference_start + 1 <= end_coord:
                            continue
                        minus[read.reference_start - start_coord + 1] += 1

                else:
                    if read.is_read1:
                        if not start_coord <= read.reference_start + 1 <= end_coord:
                            continue
                        minus[read.reference_start - start_coord + 1] += 1
                    else:
                        if not start_coord <= read.reference_end + 1 <= end_coord:
                            continue
                        plus[read.reference_end - start_coord + 1] += 1

            if libtype == "FR":
                return cls(chrom, start_coord, end_coord, plus, -minus)
            else:
                return cls(chrom, start_coord, end_coord, minus, -plus)

        except IOError:

            raise f'There is no .bam file at {bam_file_path}'

    @classmethod
    def generateobj(cls):

        return cls(None, None, None, None, None)

    def is_invalid(self):
        """
        Check the read depth object whether legal
        """
        return self.chrom is None or self.low is None or self.high is None or self.plus is None or self.minus is None

    def __add__(self, other):

        if self.is_invalid():
            return other
        if other.is_invalid():
            return self

        assert self.chrom == other.chrom, 'Cannot add depths from different chromosomes'
        assert self.low == other.low and self.high == other.high, 'Cannot add depths with different start and end ' \
                                                                  'points '
        newplus = self.plus + other.plus
        newminus = self.minus + other.minus

        return SiteDepth(self.chrom, self.low, self.high, newplus, newminus)

    def __str__(self):
        return '{0}:{1}-{2},{3},{4}'.format(self.chrom, self.low, self.high, self.plus, self.minus)
