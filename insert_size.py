#!/usr/bin/env python
###############################################################################
#
#    Print out the insert sizes of each of pair of mapped reads in a bam file
#
#    Copyright (C) 2013 Ben Woodcroft
#
###############################################################################

__copyright__ = "Copyright 2013"
__license__ = "GPLv3"
__version__ = "0.0.1"

###############################################################################

import argparse
import sys
import pysam
import os


class ReadLoader:
    """AUX: Call back for getting aligned reads

    Used in conjunction with pysam.fetch
    """
    def __init__(self):
        self.alignedReads = []

    def __call__(self, alignment):
        self.alignedReads.append(alignment)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("bam", help="BAM file to be analysed [required]")
    #parser.add_argument('positional_arg', help="Required")
    #parser.add_argument('positional_arg2', type=int, help="Integer argument")
    #parser.add_argument('positional_arg3', nargs='+', help="Multiple values")
    #parser.add_argument('--min-bases-per-read', type=int, default=300, dest='coverage_threshold', help="Amount of bases allowable per read (enriching for low-coverage contigs)")

    # parse the arguments
    args = parser.parse_args()

    pysam_object = pysam.Samfile(args.bam, 'rb')
    total_start_positions = 0
    #start_position_to_coverage_count =

    for reference, contig_length in zip(pysam_object.references, pysam_object.lengths):
        rl = ReadLoader()
        pysam_object.fetch(reference, 0, contig_length, callback = rl)
        start_positions = {}
        for read in rl.alignedReads:
            if read.is_secondary or not read.is_proper_pair:
                continue

            # Only count reads that are first in the pair and are not reverse complemented
            if read.is_read2 or read.is_reverse:
                continue

            print read.tlen



