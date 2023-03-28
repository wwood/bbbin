#!/usr/bin/env python3

###############################################################################
#
#    Copyright (C) 2023 Ben Woodcroft
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

__author__ = "Ben Woodcroft"
__copyright__ = "Copyright 2023"
__credits__ = ["Ben Woodcroft"]
__license__ = "GPL3"
__maintainer__ = "Ben Woodcroft"
__email__ = "benjwoodcroft near gmail.com"
__status__ = "Development"

import argparse
import logging
import sys
import os
import csv
from collections import OrderedDict

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    # -i
    parent_parser.add_argument('-i', '--input', help='input file', required=True)
    parent_parser.add_argument('-r', '--reference', help='reference FASTA file', required=True)

    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    # Read reference contigs FASTA
    contig_lengths = OrderedDict()
    with open(args.reference, 'r') as f:
        for line in f:
            if line.startswith('>'):
                contig = line.strip()[1:].split(" ")[0]
                contig_lengths[contig] = 0
            else:
                contig_lengths[contig] += len(line.strip())

    # Print header of VCF
    print('''##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##source=instrain
##sample=<ID=1,name=reads.1.fq.gz.bam>''')

    for contig, length in contig_lengths.items():
        print('##contig=<ID={},length={}>'.format(contig, length))

    print('''##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tinStrain''')


    # Iterate CSV - re-order references according to order in ref FASTA
    contig_to_snvs = {}
    with open(args.input, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        for row in reader:
            #             scaffold        position        position_coverage       allele_count    ref_base        con_base        var_base        ref_freq        con_freq        var_freq        A
            #        C       T       G       cryptic class
            # NZ_CP033091.2   92      7       1       A       G       A       0.0     1.0     0.0     0       0       0       7       False   SNS
            # NZ_CP033091.2   103     8       1       A       G       A       0.0     1.0     0.0     0       0       0       8       False   SNS
            # N
            
            # NZ_CP033092.2   82      .       C       T       . .       .      GT  1

            if row[0] not in contig_to_snvs:
                contig_to_snvs[row[0]] = []

            # inStrain reports 0-based values for “position”. The first base in a scaffold will be position “0”, second based position “1”, etc.

            contig_to_snvs[row[0]].append("\t".join([
                row[0], # scaffold
                str(int(row[1])+1), # position
                '.',
                row[4], # ref_base
                row[5], # con_base
                '.',
                'PASS',
                '.',
                'GT',
                '0']))

            # If there is a 2nd allele, print that too
            if int(row[3]) != 1:
                # Sometime we get 3 alleles according to the 4th column, but
                # unclear how to glean the 3rd allele from the TSV file. So just
                # ignore third, 4th etc alleles
                contig_to_snvs[row[0]].append("\t".join([
                    row[0], # scaffold
                    str(int(row[1])+1), # position
                    '.',
                    row[4], # ref_base
                    row[6], # var_base
                    '.',
                    'PASS',
                    '.',
                    'GT',
                    '0'
                ]))

    # Print SNVs
    for contig in contig_lengths.keys():
        for snv in contig_to_snvs[contig]:
            print(snv)