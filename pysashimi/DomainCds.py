#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/1/11 2:55 PM
__author__ = 'Zhou Ran'

import math
from .Constant import ANNOTATION


def bs(arr, threshold):
    """
    Binary Search

    :param arr:
    :param threshold:
    :return:
    """
    low = 0
    high = len(arr)

    while low < high:
        mid = math.floor((low + high) / 2)

        if arr[mid] == threshold:
            return mid
        elif arr[mid] < threshold and mid != low:
            low = mid
        elif arr[mid] > threshold and mid != high:
            high = mid
        else:
            high = low = low + 1

    return low


def calculateinterval(cdsinterval, txregion):
    """
    calculate the interval which were included
    :param cdsinterval:
    :param txregion:
    :return:
    """
    assert txregion[0] < txregion[1], "start site not smaller than end site in the given region"
    cdsinterval = list(map(list, sorted(cdsinterval, key=lambda x: x[0])))

    lindex, rindex = CdsDmain.domainlocation(cdsinterval, txregion)

    # remove the intersect region because one point.
    if not all([max(map(lambda x: max(x), cdsinterval)) > min(txregion),
                min(map(lambda x: max(x), cdsinterval)) < max(txregion)]):
        return None

    info = cdsinterval[lindex:rindex]


    if info[0][0] < txregion[0]:
        info[0][0] = txregion[0]

    if info[-1][-1] > txregion[1]:
        info[-1][-1] = txregion[1]

    return info


class CdsDmain:
    def __init__(self, cdsinfo, domaininfo, strand):
        """

        :param cdsinfo: a list contain the cds region information. `[(),(),()]`
        :param domaininfo: a list contain the domain information
        """
        self.cdsinfo = sorted(cdsinfo, key=lambda x: x[0])
        self.domaininfo = domaininfo
        self.strand = strand

    def __relativedistance(self):
        last = 1
        relativecds = []
        cdsinfo = self.realcdsinfo

        # if self.strand == "-":
        #     cdsinfo = self.cdsinfo[::-1]
        # print(cdsinfo)
        for num, cds_ in enumerate(cdsinfo):
            onecdsregion = self.dist(cds_)
            relativecds.append((last, last + onecdsregion))
            last += onecdsregion + 1
            # if num == 0:
            #     relativecds.append((last, last + onecdsregion))
            #     last += onecdsregion + 1
            #
            # else:
            #     relativecds.append((last, last + onecdsregion))
            #     last += onecdsregion + 1
        return relativecds

    @property
    def relativecds(self):
        """
        every length of cds interval
        :return:
        """
        return self.__relativedistance()

    @property
    def realcdsinfo(self):
        """
        cds region list
        :return:
        """
        if self.strand == "+":
            return self.cdsinfo
        return self.cdsinfo[::-1]

    @staticmethod
    def dist(lst):
        """

        :param lst: (s,e)
        :return: abs(s-e)
        """
        return abs(int(lst[0]) - int(lst[1]))

    @staticmethod
    def domainlocation(relativecds, domainregion):
        """

        :param relativecds: an element of an object of self.relativecds
        :param domainregion: (domainstart,domainend,domainname)
        :return:
        """
        r = list(map(lambda x: x[1], relativecds))
        l = list(map(lambda x: x[0], relativecds))
        lindex = bs(r, domainregion[0])
        rindex = bs(l, domainregion[1])

        return lindex, rindex

    @staticmethod
    def offsetinregion(interval, site, left=True):
        """
        here to judge the location of site, and return the offset to the first site
        :param interval: (start,end)
        :param site: site
        :return: int
        """
        l, r = interval
        # assert l <= site <= r, \
        #     "The site ({}) is not in the interval: [{}, {}]".format(site, l, r)
        if site < l:
            if left:
                return 0
            else:
                return None
        elif site > r:
            if left:
                return None
            else:
                return r - l + 1

        else:
            return site - l + 1

    @property
    def domainrelativegenomiccoordinary(self):
        """

        :return:
        """
        indexres = []
        # indexres = defaultdict()

        for d_ in self.domaininfo:
            name = d_[-1]
            lindex, rindex = self.domainlocation(self.relativecds, d_)
            # print(lindex, rindex, self.realcdsinfo)
            if lindex == rindex:
                continue
            loffset = self.offsetinregion(self.relativecds[lindex], d_[0])
            roffset = self.offsetinregion(self.relativecds[rindex - 1], d_[1], left=False)
            cdsregion = self.realcdsinfo[lindex:rindex]
            if self.strand == "+":
                u'''
                1.16 if the domain in a same region, must save a temp value
                '''
                left_orgin = list(cdsregion[0])
                right_orgin = list(cdsregion[-1])
                tmplist = list(cdsregion[0])
                tmplist[0] = left_orgin[0] + loffset
                cdsregion[0] = tuple(tmplist)

                tmplist = list(cdsregion[-1])
                tmplist[1] = right_orgin[0] + roffset
                cdsregion[-1] = tuple(tmplist)

            else:
                left_orgin = list(cdsregion[0])
                right_orgin = list(cdsregion[-1])

                tmplist = list(cdsregion[0])
                tmplist[1] = left_orgin[1] - loffset + 1
                cdsregion[0] = tuple(tmplist)

                tmplist = list(cdsregion[-1])
                tmplist[0] = right_orgin[1] - roffset
                cdsregion[-1] = tuple(tmplist)
            indexres.append((cdsregion, name))
        return indexres


def annot_classfication(domaininfo):
    """
    domaininfo was a object that like a string, for example, domainname;;domaintype
    :param domaininfo:
    :return:
    """
    domaindes, domaintype = domaininfo.split(';;')
    if domaintype == 'repeat':
        return 'repeat'
    for annotcate, detail in ANNOTATION.items():
        if domaintype in detail:
            return annotcate
    return 'other'


def main():
    # cdsinfo = [
    #     (3999560, 3999617), (4007656, 4007737), (4019070, 4019148), (4024736, 4024890), (4041888, 4042107),
    #     (4092617, 4092780), (4119668, 4119712), (4120015, 4120073), (4142612, 4142766), (4147812, 4147963),
    #     (4148612, 4148744), (4163855, 4163941), (4170205, 4170404), (4197534, 4197641), (4206660, 4206837),
    #     (4226611, 4226823), (4228443, 4228619), (4231053, 4231144), (4243133, 4243262), (4243417, 4243448),
    #     (4243543, 4243619), (4245031, 4245106), (4261527, 4261605), (4267469, 4267620), (4284766, 4284898),
    #     (4292926, 4293012), (4311270, 4311433), (4351910, 4352081), (4352202, 4352837), (4409170, 4409187)
    # ]
    cdsinfo = [(6230003, 6230073), (6233961, 6234087), (6234229, 6234311)]  # 71,127 83
    # domaininfo = [
    #     (1, 1941, '1XK-related protein 4;chain'), (334, 396, '2Helical;transmembrane region'),
    #     (424, 486, '3Helical;transmembrane region'), (733, 795, '4Helical;transmembrane region'),
    #     (907, 969, '5Helical;transmembrane region'), (982, 1044, '6Helical;transmembrane region'),
    #     (1084, 1146, '7Helical;transmembrane region'), (1177, 1245, '8Helical;transmembrane region'),
    #     (1273, 1335, '9Helical;transmembrane region'), (1360, 1422, '10Helical;transmembrane region'),
    #     (1450, 1512, '11Helical;transmembrane region'), (589, 591, '12Phosphoserine;modified residue')
    # ]
    domaininfo = [(277, 279, 'non-terminal residue;non-terminal residue')]

    a = CdsDmain(cdsinfo, domaininfo, '+')
    print(a.relativecds)
    print(a.domainrelativegenomiccoordinary)


if __name__ == '__main__':
    main()
