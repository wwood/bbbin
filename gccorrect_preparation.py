#!/usr/bin/env python
#
#    Copyright (C) 2013 Ben Woodcroft, available under GPLv3 or later
#


from optparse import OptionParser
import sys
from pprint import pprint
import pysam

class ReadLoader:
    """AUX: Call back for getting aligned reads

    Used in conjunction with pysam.fetch
    """
    def __init__(self):
        self.alignedReads = []

    def __call__(self, alignment):
        self.alignedReads.append(alignment)


if __name__ == '__main__':
    # intialise the options parser
    parser = OptionParser("\n\n %prog [options]")

    # add options hereread
    #parser.add_option("-f", "--fasta", dest="fasta", help="Fasta file of sequences to be prepped [required]")
    parser.add_option("-b", "--bam", dest="bam", help="BAM file to be analysed [required]")
    parser.add_option("-f", "--forward", dest="forward_file", help="Output forwards to this file [required]")
    parser.add_option("-r", "--reverse", dest="reverse_file", help="Output reverse reads to this file [required]")
    (opts, args) = parser.parse_args()

    sam = pysam.Samfile(opts.bam, 'rb')
    f = open(opts.forward_file,'w')
    r = open(opts.reverse_file,'w')
    for reference, contig_length in zip(sam.references, sam.lengths):
        rl = ReadLoader()
        sam.fetch(reference, 0, contig_length, callback = rl)
        print "Found",len(rl.alignedReads),"reads to consider"

        for read in rl.alignedReads:
            # Ignore unpaired reads or secondary hits - reads should only count once
            if read.is_secondary or not read.is_proper_pair:
                continue

            # Only need to work with the read1's, not their partners
            if read.is_read2:
                continue

            # Not sure how, but this appear to happen somehow. TODO: advise the user how many times
            if read.tlen < 0:
                continue

            # OK, so we have a read1. We should now be able to write out
            #TODO: check the +/- 1 are right in the reads below
            output = "1 "
            if read.is_reverse:
                output += str(read.aend-1)
            else:
                output += str(read.pos+1)

            output += " "+str(read.tlen)+"\n"

            if read.is_reverse:
                r.write(output)
            else:
                f.write(output)

    f.close()
    r.close()

    #TODO: Advise how many reads were printed out


