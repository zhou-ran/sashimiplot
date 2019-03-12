#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/1/16 3:41 PM
__author__ = 'Zhou Ran'

from collections import defaultdict


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
