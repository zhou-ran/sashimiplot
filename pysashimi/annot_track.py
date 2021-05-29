#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2021/5/29 4:25 PM

"""
here to define a class to handle all track information, including:
1. bam files
2. bigwig files
3. bed files

"""

from collections import defaultdict
import pysam
from .utils import chrom_transform

class AnnotTrack:
    """
    bed track file
    """
    def __init__(self,
                 track_dic,
                 chr,
                 tstart,
                 tend
                 ) -> None:
        # pysam.TabixFile
        self.track_dic = track_dic
        self.chr = chr
        self.tstart = int(tstart)
        self.tend = int(tend)

    @property
    def tracklist(self) -> dict:
        track_return = defaultdict(list)
        for track_label, track_file in self.track_dic.items():
            track_io = pysam.TabixFile(track_file)
            try:
                io_iter = track_io.fetch(
                    self.chr,
                    self.tstart,
                    self.tend
                    )
            except ValueError:
                io_iter = track_io.fetch(
                    chrom_transform(self.chr),
                    self.tstart,
                    self.tend
                    )
            
            for line in io_iter:
                if not line:
                    continue
                line = line.strip().split('\t')
                track_return[track_label].append(
                    [
                        int(line[1]) + 1,  # to 1-base
                        int(line[2])
                    ]
                )
                
        return track_return


def main() -> None:
    track_inf = AnnotTrack(
        {'miRNA': '/Users/zhouran/Predicted_Targets.mm10.bed.gz'},
        'chr1',
        3214499,
        3215499
    ).tracklist
    print(track_inf)


if __name__ == '__main__':
    main()
