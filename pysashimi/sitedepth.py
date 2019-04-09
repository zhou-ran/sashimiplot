#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/3/21 5:20 PM
__author__ = 'Zhou Ran'

"""
This class is for plotting the site class
"""

import pysam
import numpy as np
from .utils import checkbam
from .bamfilter import pafilter


class SiteDepth:
    def __init__(self,
                 chrm,
                 low,
                 high,
                 plus,
                 minus):
        self.low = low
        self.high = high
        self.chrm = chrm
        self.minus = minus
        self.plus = plus

    @classmethod
    def determine_depth(cls,
                        bam_file_path,
                        chrm,
                        start_coord,
                        end_coord,
                        libtype,
                        singlestrand=None,
                        readfilter=None):

        """
        calculate the site coverage based the libtype and the input strand.
        FR means: the R1 is reversed to the gene strand, and the R2 is same to the gene strand

        :param bam_file_path:
        :param chrm:
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
            relevant_reads = bam_file.fetch(reference=chrm, start=start_coord, end=end_coord)

            plus = np.zeros(end_coord - start_coord + 1, dtype='f')
            minus = np.zeros(end_coord - start_coord + 1, dtype='f')

            for read in relevant_reads:
                if readfilter:
                    if not pafilter(read, readfilter):
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

                if not start_coord <= read.reference_start + 1 <= end_coord or not start_coord <= read.reference_end + 1 <= end_coord: continue

                if read.is_reverse:

                    if read.is_read1:
                        plus[read.reference_end - start_coord + 1] += 1
                    else:
                        minus[read.reference_start - start_coord + 1] += 1

                else:
                    if read.is_read1:
                        minus[read.reference_start - start_coord + 1] += 1
                    else:
                        plus[read.reference_end - start_coord + 1] += 1

            if libtype == "FR":
                return cls(chrm, start_coord, end_coord, plus, -minus)
            else:
                return cls(chrm, start_coord, end_coord, minus, -plus)

        except IOError:

            raise 'There is no .bam file at {s}'.format(bam_file_path)

    @classmethod
    def generateobj(cls):

        return cls(None, None, None, None, None)

    def is_invalid(self):
        """
        Check the readdepth object whether legal
        """
        return self.chrm is None or self.low is None or self.high is None or self.plus is None or self.minus is None

    def __add__(self, other):

        if self.is_invalid():
            return other
        if other.is_invalid():
            return self

        assert self.chrm == other.chrm, 'Cannot add depths from different chromosomes'
        assert self.low == other.low and self.high == other.high, 'Cannot add depths with different start and end points'
        newplus = self.plus + other.plus
        newminus = self.minus + other.minus

        return SiteDepth(self.chrm, self.low, self.high, newplus, newminus)

    def __str__(self):
        return '{0}:{1}-{2},{3},{4}'.format(self.chrm, self.low, self.high, self.plus, self.minus)
