#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/1/10 4:25 PM
__author__ = 'Zhou Ran'

from collections import defaultdict
from itertools import chain
import pysam

from pyUniprot import *
from DomainCds import *


class AttrDict(dict):

    def __init__(self):
        dict.__init__(self)

    def __setattr__(self, name, value):
        self[name] = value

    def __getattr__(self, name):
        return self[name]


class GTFFeature(object):
    """
    Retrieve line from GFFfile class, and return information line by line.
    """

    def __init__(self, chrom=None, source=None, featuretype=None, start=None, end=None,
                 score=None, strand=None, phase=None, attributes=None):

        self._chrom = chrom
        self._source = source
        self._featuretype = featuretype
        self._start = start
        self._end = end
        self._score = score
        self._strand = strand
        self._phase = phase
        self._attributparse(attributes)

    def _attributparse(self, attributes):
        self.attributes = AttrDict()
        if attributes:
            items = attributes.split(';')
            for i in items:
                if len(i) == 0: continue
                try:
                    name, key = i.strip().split()
                except ValueError:
                    continue
                key = key.replace('"', '')
                setattr(self.attributes, name, key)
        return self.attributes

    @property
    def chrom(self):
        return self._chrom

    @property
    def start(self):
        return int(self._start)

    @property
    def end(self):
        return int(self._end)

    @property
    def feature(self):
        return self._featuretype

    @property
    def strand(self):
        return self._strand

    def __len__(self):
        length = int(self._end) - int(self._start) + 1
        if length < 1:
            raise ValueError('Zero- or negative length feature')
        return length

    def __repr__(self):
        return "<GTF;gene {}>".format(self.attributes["gene_id"])


class mRNA:
    def __init__(self, chr, tstart, tend, gtf, genename=None):
        """

        :param chr:
        :param estart:
        :param eend:
        :param gtf:
        """

        self.chr = chr
        self.tstart = tstart
        self.tend = tend
        self.gtf = pysam.TabixFile(gtf)
        self.genename = genename
        self._validfeature = set(("exon", "CDS"))

        self.txlst_ = self.__txlst()

    def __txlst(self):
        """

        :return: all CDS and exon were contained in the given region
        """
        resdict = defaultdict(lambda: defaultdict(list))
        lines = self.gtf.fetch(self.chr, self.tstart, self.tend)
        if not lines:
            return ''
        for line in lines:
            line = GTFFeature(*line.strip().split("\t"))
            _feature = line.feature
            line_gid = line.attributes["gene_id"]

            u'''
            support the single gene mode 
            '''

            if self.genename:
                if line_gid != self.genename:
                    continue

            if _feature in self._validfeature:
                resdict[line.attributes["transcript_id"]][_feature].append((line.start, line.end))
                resdict[line.attributes["transcript_id"]]["strand"] = line.strand
            else:
                continue
        for t, cdsexon in resdict.items():
            strand = cdsexon['strand']
            if "CDS" in cdsexon:
                domainlst = []
                uniprotinfo = Uniprot(t)
                if not uniprotinfo.domain:
                    continue
                for domain in uniprotinfo.domain:
                    # if strand == "+":
                    domainlst.append(
                        (
                            domain.start,
                            domain.end,
                            ';'.join([domain.name, domain.type])
                        )
                    )
                tmpres = CdsDmain(cdsexon["CDS"], domainlst, strand).domainrelativegenomiccoordinary
                resdict[t]["domain"] = tmpres
        return resdict

    # [(1, 1257, 'Transcription factor SOX-17;chain'), (832, 1254, 'Sox C-terminal;domain'),
    #  (202, 408, 'HMG box;DNA-binding region'), (934, 1038, 'Gln/Pro-rich;compositionally biased region'),
    #  (964, 978, 'Poly-Pro;compositionally biased region'), (1, 384,  'In isoform 2.;splice variant'),
    #  (220, 267, 'helix;helix'), (271, 279, 'strand;strand'), (283, 324, 'helix;helix'), (331, 393, 'helix;helix')]

    @property
    def txlst(self):
        txlist = []
        if not self.txlst_:
            return ''
        for t, cdsexon in self.txlst_.items():
            txlist.append({t: cdsexon})
        return txlist

    @property
    def maxmin(self):
        """
        collapse the nest dict
        :return:
        """
        if not self.txlst_:
            return ''

        d = list(chain(*list(map(lambda x: list(x.values())[0], self.txlst_.values()))))
        return d

    @staticmethod
    def _returnregion(regiondict, type_):
        """

        :return: all region for a single type annotation
        """

        mRNAlist = []
        for i in regiondict:
            for k, v in i.items():
                mRNAlist.append(v[type_])
        return mRNAlist

    @property
    def exon(self):
        """
        :return: [[(s1,e1)],[(s1,e1)]]
        """
        return self._returnregion(self.txlst, 'exon')

    @property
    def exonstarts(self):
        return list(map(lambda x: x[0], chain(*self.exon)))

    @property
    def exonend(self):
        return list(map(lambda x: x[1], chain(*self.exon)))

    @property
    def domain(self):
        """
        :return: [[(s1,e1)],[(s1,e1)]]
        """
        return self._returnregion(self.txlst, 'domain')

    @property
    def cds(self):
        """
        :return: [[(s1,e1)],[(s1,e1)]]
        """
        return self._returnregion(self.txlst, 'cds')

    @staticmethod
    def domainsite(tdic, param):
        """
        :param tdic: a dict contain the CDS,exon and strand information
        :return:
        """
        if tdic['strand'] == "+":
            end_ = max(map(lambda x: x[0], tdic[param]))
            start_ = min(map(lambda x: x[0], tdic[param]))
            return end_, start_

        start_ = max(map(lambda x: x[1], tdic[param]))
        end_ = min(map(lambda x: x[1], tdic[param]))
        return end_, start_

    def max_min(self):
        """
        here to return the all start (end) site in CDS and exon
        :return: (minsite,maxsite)
        """
        d = map(lambda x: list(x.values())[0], self.txlst_.values())
        d = list(chain(*d))
        left, right = zip(*d)
        return left, right


def main(file):
    tre = mRNA('18', 20684599, 20746404, file, genename='ENSMUSG00000056124')
    print(tre.exonstarts)

    mRNAlist = []
    for i in tre.txlst:
        for k, v in i.items():
            mRNAlist.append(v['exon'])


if __name__ == '__main__':
    import logging
    import sys

    logger = logging.getLogger()
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        '%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG)

    main(sys.argv[1])
