#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/4/8 10:58 PM
__author__ = 'Zhou Ran'

"""
This code for load the machine learning model, and add the probability distribution for the pa site in the cluster
"""

import re

import numpy as np
import pysam


def DNA_matrix(seq):
    tem2 = ['[aA]', '[cC]', '[gG]', '[tT]']
    for i in range(len(tem2)):
        ind = [m.start() for m in re.finditer(tem2[i], seq)]
        tem = np.zeros(len(seq), dtype=np.int)
        tem[ind] = 1
        if i == 0:
            a = np.zeros((len(seq), 4))
        a[..., i] = tem
    return a


def DNA_matrix_(seq, exprs):
    tem2 = ['[aA]', '[cC]', '[gG]', '[tT]']
    for i in range(len(tem2)):
        ind = [m.start() for m in re.finditer(tem2[i], seq)]
        tem = np.zeros(len(seq), dtype=np.int)
        tem[ind] = 1
        if i == 0:
            a = np.zeros((len(seq) + 1, 4))
        a[..., i] = np.append(tem, exprs)
    return a


def revseq(seq):
    seqdic = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
        "N": "N"
    }
    res = ''.join([seqdic[i] for i in list(seq)])
    return res


def fetchseq(chrom, s, fa, strand):
    if strand == '+':
        seq = fa.fetch(chrom, s - 40, s + 30)
    else:
        seq = fa.fetch(chrom, s - 30, s + 40)

    if strand == '-':
        seq = revseq(seq)[::-1]
    return seq


def modelrate_(chrom, peak_s, peak_e, plotregion_s, plotregion_e, strand, model, genemtx=None):
    from keras.models import load_model
    fa = pysam.FastaFile('/mnt/raid61/Microwell/mm10/fasta/genome.fa')
    seqs = []
    for i in range(int(peak_s), int(peak_e)):
        seqs.append(fetchseq(chrom, i, fa, strand))

    for i in range(len(seqs)):
        tem = seqs[i].rstrip()
        if i == 0:
            dat_x = np.zeros((len(seqs), len(tem) + 1, 4))
        dat_x[i,] = DNA_matrix_(tem, genemtx[i])

    prob = np.zeros(plotregion_e - plotregion_s + 1, dtype='f')

    classifier = load_model(model)

    pred_ = classifier.predict(dat_x)

    # print(pred_[:, 1].shape)
    # print(prob.shape)
    # print(prob[peak_s - plotregion_s + 1:peak_e - plotregion_s + 1])
    # print(pred_[:, 1])
    prob[peak_s - plotregion_s:peak_e - plotregion_s] = pred_[:, 1]

    return prob


#
#
def modelrate(chrom, peak_s, peak_e, plotregion_s, plotregion_e, strand, model):
    from keras.models import load_model
    fa = pysam.FastaFile('/mnt/raid61/Microwell/mm10/fasta/genome.fa')
    seqs = []
    for i in range(int(peak_s), int(peak_e)):
        seqs.append(fetchseq(chrom, i, fa, strand))

    for i in range(len(seqs)):
        tem = seqs[i].rstrip()
        if i == 0:
            dat_x = np.zeros((len(seqs), len(tem), 4))
        dat_x[i,] = DNA_matrix(tem)

    prob = np.zeros(plotregion_e - plotregion_s + 1, dtype='f')

    classifier = load_model(model)

    pred_ = classifier.predict(dat_x)

    # print(pred_[:, 1].shape)
    # print(prob.shape)
    # print(prob[peak_s - plotregion_s + 1:peak_e - plotregion_s + 1])
    # print(pred_[:, 1])
    prob[peak_s - plotregion_s:peak_e - plotregion_s] = pred_[:, 1]

    return prob
