#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
from collections import defaultdict

import numpy
import pysam

from .cigarutils import fetch_intron
from .bamfilter import pafilter
from .utils import checkbam

logger = logging.getLogger("MAIN")

"""
This script were migrated from spliceplot
"""


class ReadDepth:
    def __init__(self,
                 chrm,
                 low,
                 high,
                 wiggle,
                 junctions_dict):

        self.low = low
        self.high = high
        self.chrm = chrm
        self.wiggle = wiggle
        self.junctions_dict = junctions_dict

    @classmethod
    def determine_depth(cls,
                        bam_file_path,
                        chrm,
                        start_coord,
                        end_coord,
                        readFilter=None):
        """
        calculate the coverage at each base between start_coord and endcoord.

        :param bam_file_path:
        :param chrm:
        :param start_coord: int
        :param end_coord: int
        :return: Numpy array
        """

        # checkbam(bam_file_path)
        try:
            bam_file = pysam.Samfile(bam_file_path, 'rb')
            relevant_reads = bam_file.fetch(reference=chrm, start=start_coord, end=end_coord)

            depth_vector = numpy.zeros(end_coord - start_coord + 1, dtype='f')
            spanned_junctions = defaultdict(int)

            for read in relevant_reads:
                if readFilter:
                    if not pafilter(read, readFilter):
                        continue

                # make sure that the read can be used
                cigar_string = read.cigar

                # each read must have a cigar string
                if cigar_string == None:
                    continue
                #
                # # read cannot have insertions or deletions
                # contains_indel = False

                # for cigar_event in cigar_string:
                #     if cigar_event[0] == 1 or cigar_event[0] == 2:
                #         contains_indel = True
                #         break
                #
                # if contains_indel:
                #     continue

                for index, base_position in enumerate(read.positions):
                    base_position += 1
                    if start_coord <= base_position <= end_coord:
                        depth_vector[base_position - start_coord] += 1

                intronbound = fetch_intron(read.reference_start, cigar_string)
                if intronbound:
                    for intronbound_ in intronbound:
                        junction_name = '{}:{}-{}'.format(chrm,
                                                          intronbound_[0],
                                                          intronbound_[1]
                                                          )
                        spanned_junctions[junction_name] += 1

            return cls(chrm, start_coord, end_coord, depth_vector, spanned_junctions)
        except IOError:

            raise 'There is no .bam file at {0}'.format(bam_file_path)

    @classmethod
    def generateobj(cls):

        return cls(None, None, None, None, None)

    def is_invalid(self):
        """
        Check the readdepth object whether legal
        """
        return self.chrm is None or self.low is None or self.high is None or self.wiggle is None or self.junctions_dict is None

    def __add__(self, other):
        """
        if self or other is empty, then return the other
        :param other:
        :return:
        """

        if self.is_invalid():
            return other
        if other.is_invalid():
            return self

        assert self.chrm == other.chrm, 'Cannot add depths from different chromosomes'
        assert self.low == other.low and self.high == other.high, 'Cannot add depths with different start and end points'
        new_wiggle = self.wiggle + other.wiggle

        new_junctions_dict = defaultdict(int)

        # for key, value in self.junctions_dict.items():
        #     if key in other.junctions_dict:
        #         new_junctions_dict[key] = value + other.junctions_dict[key]
        #     else:
        #         new_junctions_dict[key] = value
        u''' 2019.3.27 not intersect, but union
        '''

        for junc_key in set(list(self.junctions_dict.keys()) + list(other.junctions_dict.keys())):
            try:
                new_junctions_dict[junc_key] += self.junctions_dict[junc_key]
            except KeyError:
                pass

            try:
                new_junctions_dict[junc_key] += other.junctions_dict[junc_key]
            except KeyError:
                pass

        return ReadDepth(self.chrm, self.low, self.high, new_wiggle, new_junctions_dict)

    def __str__(self):
        return '{0}:{1}-{2},{3},{4}'.format(self.chrm, self.low, self.high, self.wiggle, self.junctions_dict)
