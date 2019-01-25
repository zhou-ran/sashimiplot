#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/1/9 10:48 AM
__author__ = 'Zhou Ran'
"""
This is the mygene.info's wrapper to process every things that you want to do
"""

from biothings_client import get_client
import logging

logger = logging.getLogger("MAIN")


# gene_client.query('ensembl.transcript:ENSMUST00000208660', fields='symbol,name,uniprot')
# gene_client = get_client('gene')

class ClientAttributionError(Exception):
    pass


class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


class Myinfo:
    def __init__(self, iquery, fields, annotfunc):
        """

        :param iquery: query input information
        :param ofields: annotation output
        :param annotfunc: function
        """
        # 在这里添加判断iquery是list还是str，进而进行判断是否是多个list的请求

        self._variant = set(['gene',
                             'variant',
                             'chem',
                             'disease',
                             'taxon'])

        self.annotfunc = annotfunc
        self.iquery = iquery
        self.ofields = fields

        if self.annotfunc not in self._variant:
            raise ClientAttributionError(
                "The function \"{}\" can't find in \"{}\"".format(self.annotfunc, '; '.join(self._variant)))

        self._myclient = get_client(self.annotfunc)
        self.query = self.__query()

    def __query(self):
        res = self._myclient.query(self.iquery, fields=self.ofields)['hits']
        return res

    @property
    def uniquequery(self):

        u'''
         stop the process, and raise a contributio error for these un-unique hits.
        '''

        assert len(self.query) == 1, "The query information was not unique, pls check {}".format(self.iquery)
        return self.query[0]

    @property
    def uniprotid(self):
        """

        :return: uniprotID
        """
        try:
            return self.uniquequery['uniprot']['Swiss-Prot']
        except:
            return ''

    @property
    def ensembl(self):
        """

        :return: all ensemble information
        """
        try:
            return AttrDict(self.uniquequery['ensembl'])
        except:
            logging.warn("No ensembl information found!")
            return ''

    @property
    def loc(self):
        """

        :return: gene coordinary iformation,
        {
        'chr': '1',
        'end': 4409241,
        'ensemblgene': 'ENSMUSG00000025900',
        'start': 3999557,
        'strand': -1
        }
        """
        try:
            return AttrDict(self.uniquequery['genomic_pos'])
        except ValueError:
            return AttrDict(self.uniquequery['genomic_pos'][0])


def main():
    a = Myinfo('ensembl.gene:ENSG00000134905', 'all', 'gene')
    # loc = a.loc
    print(a.loc)


if __name__ == '__main__':
    main()
