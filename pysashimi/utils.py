#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/1/16 3:41 PM
__author__ = 'Zhou Ran'


def readbamlist(bamlists):
    """

    Read the bam config files
    :param bamlists: file
    :return: a dict contain all sample label and bam file path
    """

    bamdict = {}
    try:
        with open(bamlists) as fi:
            for line in fi.readlines():
                label, bampath = line.strip().split('\t')
                assert label not in bamdict, "Repeated name \'{}\' were found in your config file".format(bamlists)
                bamdict[label] = bampath
        return bamdict
    except IOError:
        raise "{} file can't be found".format(bamlists)
