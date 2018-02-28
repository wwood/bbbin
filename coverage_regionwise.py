#!/usr/bin/env python

__author__ = "Ben Woodcroft"
__copyright__ = "Copyright 2017"
__credits__ = ["Ben Woodcroft"]
__license__ = "GPL3+"
__maintainer__ = "Ben Woodcroft"
__email__ = "b.woodcroft near uq.edu.au"
__status__ = "Development"

import argparse
import logging
import sys
import os
import re

import extern


if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser()
    parent_parser.add_argument('--debug', help='output debug information',
                               action="store_true")
    parent_parser.add_argument('--quiet', help='only output errors',
                               action="store_true")

    parent_parser.add_argument('--bam_file', help='sorted bam files',
                               required=True)
    parent_parser.add_argument('--interval', type=int,
                               help='report coverage for each region this size',
                               required=True)
    args = parent_parser.parse_args()


    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel,
                        format='%(asctime)s %(levelname)s: %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p')

    # Read headers to get lengths of each contig
    headers = extern.run("samtools view -H '%s'" % args.bam_file)
    contig_to_length = {}
    for line in headers.splitlines():
        splits = line.split("\t")
        if splits[0] == '@SQ':
            if len(splits) != 3:
                raise Exception("Unexpected sequence header %s" % line)
            ref = re.sub('^SN:', '', splits[1])
            if ref in contig_to_length:
                raise Exception("Contig names cannot be duplicated, found one %s" % ref)
            contig_to_length[ref] = int(re.sub('^LN:','',splits[2]))

    cmd = "samtools view -f2 -F3852 '%s'" % args.bam_file
    out = extern.run(cmd) #TODO: stream the stdout

    def finish_a_contig(contig_lengths, refname, last_position, interval, last_count):
        length = contig_lengths[refname]
        while last_position+interval <= length:
            print "\t".join([refname, str(last_position), str(last_position+interval-1), str(last_count)])
            last_count = 0
            last_position += interval

    last_ref = None
    last_position = None
    current_total_pairs = 0
    interval = args.interval
    for line in out.splitlines():
        splits = line.split("\t") #TODO: use csv for faster
        logging.debug("Interrogating line %s" % str(splits))
        if len(splits) <= 11: raise Exception("unexpected number of fields in sam line %s" % line)
        ref = splits[2]
        start = int(splits[3])
        tlen = int(splits[8])
        if tlen < 0: continue

        # if this is a new contig
        if last_ref != ref:
            # print the last contig
            if last_ref is not None:
                finish_a_contig(contig_to_length, last_ref, last_position, interval, current_total_pairs)
            # initialise the current contig
            current_total_pairs = 0
            last_position = 1
            last_ref = ref

        if last_position > start:
            raise Exception("BAM file not sorted?")
        # if read pair is completely contained within the current interval, increment the count
        if last_position+interval > start:
            if last_position+interval >= start+tlen:
                logging.debug("Taking it")
                current_total_pairs += 1
            else:
                logging.debug("skipping it")
                # else if start is in the current interval, ignore this pair
                pass
        else:
            logging.debug("Going to next -------------------------------------------------")
            # else the start is in the next interval, if not further
            # print the last known interval
            # print the empty intervals in between
            while last_position+interval <= start:
                print "\t".join([last_ref, str(last_position), str(last_position+interval-1), str(current_total_pairs)])
                current_total_pairs = 0
                last_position += interval
            current_total_pairs = 0
            # initialise the current interval
            if last_position+interval > start+tlen:
                # read pair completely contained
                current_total_pairs = 1

    # print the last contig
    if last_ref is not None:
        finish_a_contig(contig_to_length, last_ref, last_position, interval, current_total_pairs)
