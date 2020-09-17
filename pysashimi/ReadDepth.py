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
                 chrom,
                 low,
                 high,
                 wiggle,
                 junctions_dict):

        self.low = low
        self.high = high
        self.chrom = chrom
        self.wiggle = wiggle
        self.junctions_dict = junctions_dict

    @classmethod
    def determine_depth(cls,
                        bam_file_path,
                        chrom,
                        start_coord,
                        end_coord,
                        scale=False,
                        barcode_info=None,
                        cell_tag=None,
                        umi_tag=None,
                        readFilter=None):
        """
        calculate the coverage at each base between start_coord and endcoord.

        :param bam_file_path:
        :param chrom:
        :param start_coord:
        :param end_coord:
        :param scale:
        :param readFilter:
        :return:
        """

        try:
            bam_file = pysam.Samfile(bam_file_path, 'rb')

            relevant_reads = bam_file.fetch(reference=chrom, start=start_coord, end=end_coord)

            depth_vector = numpy.zeros(end_coord - start_coord + 1, dtype='f')
            spanned_junctions = defaultdict(int)

            if barcode_info:
                cluster_cov = defaultdict(lambda :defaultdict())
                for cluster in set(barcode_info.values()):
                    cluster_cov[cluster]['depth'] = depth_vector
                    cluster_cov[cluster]['junc'] = spanned_junctions

            for read in relevant_reads:
                if barcode_info:
                    try:
                        current_bc = read.get_tag(cell_tag)
                        _ = read.get_tag(umi_tag)
                        current_bc = current_bc.strip().split('-')[0]
                        if current_bc not in barcode_info:
                            continue
                        else:
                            current_cluster = barcode_info[current_bc]
                    except KeyError:
                        continue

                if readFilter:
                    if not pafilter(read, readFilter):
                        continue

                # make sure that the read can be used
                cigar_string = read.cigar

                # each read must have a cigar string
                if cigar_string is None:
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
                        if barcode_info:
                            cluster_cov[current_cluster]['depth'][base_position - start_coord] += 1
                        else:
                            depth_vector[base_position - start_coord] += 1

                intronbound = fetch_intron(read.reference_start, cigar_string)
                if intronbound:
                    for intronbound_ in intronbound:
                        junction_name = '{}:{}-{}'.format(chrom,
                                                          intronbound_[0],
                                                          intronbound_[1]
                                                          )
                        if barcode_info:
                            cluster_cov[current_cluster]['junc'][junction_name] += 1
                        else:
                            spanned_junctions[junction_name] += 1
            if scale and not barcode_info:
                max_value = numpy.max(depth_vector)
                if max_value != 0:
                    depth_vector = depth_vector/max_value * 100
            if barcode_info:
                cluster_cov_return = {}
                for cluster, dep_junc_info in cluster_cov.items():
                    cluster_cov_return[cluster] = cls(
                        chrom,
                        start_coord,
                        end_coord,
                        depth_vector['depth'],
                        depth_vector['junc']
                    )
                return cluster_cov_return

            return cls(chrom, start_coord, end_coord, depth_vector, spanned_junctions)

        except IOError:

            raise Exception(f'There is no .bam file at {bam_file_path}')


    @classmethod
    def generateobj(cls):

        return cls(None, None, None, None, None)

    def is_invalid(self):
        """
        Check the readdepth object whether legal
        """
        return self.chrom is None or self.low is None or self.high is None or self.wiggle is None or self.junctions_dict is None

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

        assert self.chrom == other.chrom, 'Cannot add depths from different chromosomes'
        assert self.low == other.low and self.high == other.high, \
            'Cannot add depths with different start and end points'
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

        return ReadDepth(self.chrom, self.low, self.high, new_wiggle, new_junctions_dict)

    def __str__(self):
        return '{0}:{1}-{2},{3},{4}'.format(self.chrom, self.low, self.high, self.wiggle, self.junctions_dict)
