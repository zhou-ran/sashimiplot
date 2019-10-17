#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/5/8 9:08 PM
__author__ = 'Zhou Ran'

import numpy as np
from loguru import logger
from scipy import stats


class EM:
    """
    EM algorithm for analyzing the pA site information
    """

    def __init__(self,
                 siteCounts,
                 alpha=3,
                 peakStrand=None,
                 iterTime=500):
        self.siteCounts = np.array(siteCounts)
        self.peakStrand = peakStrand
        self.iterTime = iterTime
        self.alpha = alpha

        self.__peakLength = len(self.siteCounts)
        self.__threeQuarter = np.percentile(self.siteCounts, 75)

        self._index = self.generateIndex()
        self._mtx = self.generateArray()
        self.pi = self.initPi()

    @staticmethod
    def pDensity(x, _lambda):
        """
        the exp distribution for simulating the distribution of all reads for a real peak.
        :param x:
        :param _lambda:
        :return:
        """
        return stats.expon.pdf(x,
                               scale=1 / _lambda
                               )

    def generateIndex(self):
        """
        generate the location index for calculating the offset of distance to the pA site.
        :param length:
        :return:
        """

        return np.linspace(1,
                           self.__peakLength,
                           num=self.__peakLength
                           )

    def generateArray(self):
        """
        generate the matrix in E step
        :param length: int, means there was `length` sites
        :return:
        """

        return np.zeros([
            self.__peakLength,
            self.__peakLength
        ]
        )

    def initPi(self):
        """
        init the pi, which was the prior probability for EM algorithm
        :return:
        """

        return np.ones(self.__peakLength) / self.__peakLength

    def loglikelihood(self):
        return np.sum(np.log(np.sum(self._mtx, axis=0) * self.pi + 1))

    def eStep(self):
        for _index in range(self.__peakLength):
            distanceOffset = _index - self._index + 1
            self._mtx[_index, :] = self.pi[_index] * self.pDensity(distanceOffset, 1 / 5)

    def mStep(self):
        self._mtx = self._mtx / np.sum(self._mtx, axis=0)

        self._mtx = np.nan_to_num(self._mtx)

        tmpMtx = self.siteCounts * self._mtx
        tmpMtx = np.apply_along_axis(lambda x: max(0, sum(x) - self.__threeQuarter), 1, tmpMtx)

        self.pi = tmpMtx / np.sum(tmpMtx)

    def fit(self, tol=1e-4):
        c_lllh = 0
        for i in range(self.iterTime):
            self.eStep()
            # if i >= 1:
            #     self.alpha = self.alpha * (i - 1) /
            lllh = self.loglikelihood()
            self.mStep()

            if c_lllh == 0:
                c_lllh = lllh
            else:
                if lllh - c_lllh < tol:
                    logger.info('Convergence after {} cycling'.format(i))
                    break
                else:
                    c_lllh = lllh
