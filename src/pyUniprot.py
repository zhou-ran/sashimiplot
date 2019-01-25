#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/1/10 2:04 PM
__author__ = 'Zhou Ran'
import logging
from collections import namedtuple
import requests as rq
import xmltodict

logger = logging.getLogger("MAIN")


class ErrorCrawl(Exception):
    pass


# TODO,add a function, if not a protein_coding, continue
class Uniprot:

    def __init__(self, uniprot_id, frmt="XML", database="uniprot"):
        self.ui = uniprot_id
        self.frmt = frmt
        self.database = database

        self._validfrmt = ['txt', 'XML', 'rdf', 'gff', 'fasta']
        self._url = 'https://www.uniprot.org'
        self.info = self.__info()

    def __requesturl(self):
        """
        request the uniprot url and return a xml object
        :return:
        """
        assert self.frmt in set(self._validfrmt), \
            'The {} \'format\' is not legal, must be one of {}'.format(self.frmt, ','.join(self._validfrmt))

        url = "{}/uniprot/?query={}&format={}".format(self._url, self.ui, self.frmt)

        try:
            result = rq.get(url, timeout=10)

        except ConnectionError:
            raise "Can't crawl {}, pls check your network".format(url)

        return result

    def __info(self):
        xmldic = xmltodict.parse(self.__requesturl().text, attr_prefix='', cdata_key='')['uniprot']['entry']

        return xmldic

    @property
    def feature(self):
        try:
            f = self.info['feature']
            res = []
            if not isinstance(f, list):
                res.append(f)
                return res
            return self.info['feature']

        except:
            return None

    @property
    def domain(self):
        if not self.feature:
            return None
        domaininfo = namedtuple('domaininfo', ['name', 'type', 'start', 'end'])
        domainres = []
        for d_ in self.feature:
            if 'description' in d_:
                description = 'description'
                type_ = 'type'
            else:
                description = 'type'
                type_ = 'type'

            if 'position' in d_['location']:
                domainres.append(domaininfo._make([d_[description],
                                                   d_[type_],
                                                   (int(d_['location']['position']['position']) - 1) * 3 + 1,
                                                   int(d_['location']['position']['position']) * 3]))
            else:
                try:
                    domainres.append(domaininfo._make([d_[description],
                                                       d_[type_],
                                                       (int(d_['location']['begin']['position']) - 1) * 3 + 1,
                                                       int(d_['location']['end']['position']) * 3]))
                except KeyError:
                    '''
                    1.24 there was a location error, because the location status was unkown!!
                    '''
                    if 'status' in d_['location']['begin'] or 'status' in d_['location']['end']:
                        continue

        return domainres

    @property
    def dbreference(self):
        """
        OrderedDict([('type', 'PROSITE'),
        ('id', 'PS50095'), ('property',
        [OrderedDict([('type', 'entry name'),('value', 'PLAT')]),
        OrderedDict([('type', 'match status'), ('value', '6')])])])
        :return: dict in list
        """
        return self.info['dbReference']

    @property
    def ensemblinfo(self):
        """
        OrderedDict([('type', 'Ensembl'), ('id', 'ENSMUST00000208660'),
        ('property', [OrderedDict([('type', 'protein sequence ID'),
        ('value', 'ENSMUSP00000146439')]),
        OrderedDict([('type', 'gene ID'), ('value', 'ENSMUSG00000025900')])])])
        :return: ensemblinfo
        """
        for info in self.dbreference:
            if info['type'] == 'Ensembl':
                return info['id']

        return ""


def main(id):
    print('main')
    b = Uniprot(id)
    print(b.feature)
    print('done')
    print(b.ensemblinfo, b.domain)


if __name__ == '__main__':
    import sys

    main(sys.argv[1])
