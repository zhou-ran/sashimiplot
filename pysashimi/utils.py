#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/1/16 3:41 PM
__author__ = 'Zhou Ran'
import os
import logging
from collections import defaultdict
import pysam

def readbamlist(bamlists):
    """

    Read the bam config files
    :param bamlists: file
    :return: a dict contain all sample label and bam file path
    """

    bamdict = defaultdict(list)
    try:
        with open(bamlists) as fi:
            for line in fi.readlines():
                try:
                    label, bampath = line.strip().split('\t')
                    bamdict[label].append(bampath)

                except ValueError:
                    line = line.strip().split('\t')
                    bamdict[line[0]].extend(line[1:])

        return bamdict
    except IOError:
        raise "{} file can't be found".format(bamlists)


def checkbam(bamfile):
    """
    check bam index
    :param file:
    :return:
    """
    if not os.path.exists(bamfile + '.bai'):
        logging.info("Index the bam file with 4 cores")
        pysam.index(bamfile, "-@", "4")
    else:
        # Check the index file whether is elder than bam file, if not re-generate
        if os.path.getmtime(bamfile) > os.path.getmtime(bamfile + '.bai'):
            logging.info('The index file is older than %s file, removing the index file' % bamfile)
            os.remove(bamfile + '.bai')
            logging.info("Index the bam file with 4 cores")
            pysam.index(bamfile, "-@", "4")
        else:
            logging.info("Index of the bam file was complete!")
