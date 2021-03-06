#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/1/16 3:41 PM
__author__ = 'Zhou Ran'
import os
import sys
from loguru import logger
from collections import defaultdict
import pysam


def readbamlist(bamlists):
    """

    Read the bam config files
    :param bamlists: file
    :return: a dict contain all sample label and bam file path
    """

    bamdict = defaultdict(list)
    colordict = defaultdict(set)
    try:
        with open(bamlists) as fi:
            for line in fi.readlines():
                try:
                    line = line.strip().split('\t')
                    if len(line) == 2:
                        label, bampath = line
                        bamdict[label].append(bampath)
                    else:
                        label, bampath, color = line
                        colordict[label].add(color)
                        bamdict[label].append(bampath)
                except ValueError:
                    raise "There was more than three or only one columns in {} file".format(bamlists)

        return bamdict, colordict
    except IOError:
        raise "{} file can't be found".format(bamlists)


def checkbam(bamfile):
    """
    check bam index
    :param file:
    :return:
    """
    if not os.path.exists(bamfile + '.bai'):
        logger.info("Index the bam file with 4 cores")
        pysam.index(bamfile, "-@", "4")
    else:
        # Check the index file whether is elder than bam file, if not re-generate
        if os.path.getmtime(bamfile) > os.path.getmtime(bamfile + '.bai'):
            logger.info('The index file is older than %s file, removing the index file' % bamfile)
            try:
                os.remove(bamfile + '.bai')
                logger.info("Index the bam file with 4 cores")
                pysam.index(bamfile, "-@", "4")
            except PermissionError:
                logger.info("Index time error, but no permission to process, Skip")

        else:
            logger.info("Index of the bam file was complete!")


def read_track(trackfile):
    """

    Read the track config files
    :param bamlists: file
    :return: a dict contain all sample label and target tracking information
    """

    trackdict = defaultdict()
    try:
        with open(trackfile) as fi:
            for line in fi.readlines():
                line = line.strip().split('\t')
                if len(line) == 2:
                    label, track_path = line
                    trackdict[label] = track_path
                else:
                    sys.exit(
                        "There was more than more or only one columns in {} file".format(trackfile)
                        )
        return trackdict
    except IOError:
        raise "{} file can't be found".format(trackfile)


def chrom_transform(chromosome_id):
    """
    convert chromosome id from with or without `chr` suffix
    """
    if chromosome_id.startswith('chr'):
        if chromosome_id == 'chrM':
            chrom_id = 'chrM'
        else:
            chrom_id = chromosome_id.replace('chr', '')
    else:
        if chromosome_id == 'chrM':
            chrom_id = 'MT'
        else:
            chrom_id = f'chr{chromosome_id}'

    return chrom_id