#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/1/17 3:27 PM
__author__ = 'Zhou Ran'

COLOR = ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
         "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
         "#E41A1C"]

# COLOR = [
#     '#8dd593',
#     '#c6dec7',
#     '#ead3c6',
#     '#f0b98d',
#     '#ef9708',
#     '#0fcfc0',
#     '#9cded6',
#     '#d5eae7',
#     '#f3e1eb',
#     '#f6c4e1',
#     '#f79cd4'
# ]

# This dataset were download from intronic pA destipted in nature communication.

DOMAINFILTER = {"active site", "domain", "signal peptide", "transmembrane region", "repeat", "zinc ﬁnger region",
                "compositionally biased region", "DNA-binding region", "region of interest",
                "lipid moiety-binding region", "short sequence motif", "calcium-binding region",
                "nucleotide phosphate-binding region", "metal ion-binding site", "topological domain"}

DNABIND = {'C2H2-type', 'PHD-type', 'C3H1-type', 'KRAB', 'Bromo', 'Chromo', 'DNA-binding', 'C4-type', 'CHCR',
           'A.T hook', 'bZIP', 'bHLH', 'CCHC-type', 'CHCH', 'Bromodomain-like', 'CH1', 'C6-type', 'A.T hook-like',
           'C4H2 - type', 'CHHC-type'}

ACTIVE = {'active site', 'catalytic sites'}

TRANSREGION = {'transmembrane region', 'ABC transmembrane type-1', 'ABC transporter', 'ABC transmembrane type-2'}

PPI = {"WD", "ANK", "TPR", "LRR", "HEAT", "Sushi", "EF-hand", "ARM", "PDZ", "PH", "SH3", "RING-type",
       "LIM zinc-binding", "WW", "SH2", "BTB", "FERM", "CH", "Rod", "Coil 1A", "MH2", "WD40-like repeat",
       "t-SNARE coiled-coil homology", "Coil 1B", "Cbl-PTB", "Coil", "CARD", "SH2-like", "DED", "IRS-type PTB",
       "SP-RING-type", "EF-hand-like", "RING-CHtype", "v-SNARE coiled-coil homology", "Arm domain",
       "LIM protein-binding", "GYF", "PDZ domain-binding", "and PDZD11-binding"}

RNABIND = {"RRM", "SAM", "KH", "DRBM", "RBD", "Piwi", "PAZ", "S1 motif", "Pumilio", "THUMP"}

ANNOTATION = {
    "PPI": PPI,
    "TMD": TRANSREGION,
    "activesite": ACTIVE,
    "dnabinding": DNABIND,
    "rnabinding": RNABIND
}

