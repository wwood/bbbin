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
import polars as pl

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser()
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--input-taxonomy', help='2 column <ID><tab><taxonomy> with 7 level taxonomy', required=True)
    # This idea currently buggy
    parent_parser.add_argument('--no-add-prefixes', help='Do not add prefixes to taxons', action="store_true")

    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')

    # Read taxonomy map
    df = pl.read_csv(args.input_taxonomy, separator='\t', has_header=False)
    df.columns = ['id','taxonomy']

    levels = ['d__','p__','c__','o__','f__','g__','s__']

    # Gather known taxonomies. This checks for duplicates as well as providing a
    # data structure to impute dummy middle rank taxons.
    taxon_to_parent = {}
    for row in df.rows(named=True):
        id = row['id']
        tax = row['taxonomy']
        taxons = list([t.strip().replace(' ','_') for t in tax.split(';')])
        if len(levels) != len(taxons):
            raise Exception("Unexpected taxon format: %s" % row['taxonomy'])

        for i in range(1, len(taxons)):
            if args.no_add_prefixes: # This idea currently buggy
                taxon = taxons[i]
                parent = taxons[i-1]
            else:
                taxon = levels[i]+taxons[i]
                parent = levels[i-1]+taxons[i-1]
            if len(taxon) <= 3: # ignore this missing rank, this will be filled by lower ranks
                if i == 6:
                    raise Exception("Taxon %s has no species" % tax)
                continue
            elif len(parent) <= 3:
                # Add a dummy entry for this parent, plus all other parents
                # which are missing above that
                j = i-1
                last_taxon = taxon
                while j >= 0 and taxons[j] == '':
                    parent_name = levels[j] + 'FILLED_' + taxon
                    if last_taxon in taxon_to_parent and taxon_to_parent[last_taxon] != parent_name:
                        raise Exception("Taxon %s has multiple parents: %s and %s" % (last_taxon, taxon_to_parent[last_taxon], parent_name))
                    taxon_to_parent[last_taxon] = parent_name
                    last_taxon = parent_name
                    j -= 1
                if args.no_add_prefixes:
                    taxon_to_parent[last_taxon] = taxons[j]
                else:
                    taxon_to_parent[last_taxon] = levels[j] + taxons[j]
            elif taxon not in taxon_to_parent:
                taxon_to_parent[taxon] = parent
            elif taxon_to_parent[taxon] != parent:
                raise Exception("Taxon %s has multiple parents: %s and %s" % (taxon, taxon_to_parent[taxon], parent))
            # Otherwise, taxon already has the correct parent, nothing new or
            # untoward here

    for row in df.rows(named=True):
        # matching rename_proteomes.py
        # genome = os.path.basename(proteome.split('.')[0])
        # # Use species names which do not contain underscores, so no speciesRax mapping file is required
        # genome_name = genome.replace('_','')
        id2 = row['id'].split('.')[0].replace('_','')
        # Build the taxonomy from the species up
        current_taxon = row['taxonomy'].split(';')[-1].strip().replace(' ','_')
        if not args.no_add_prefixes:
            current_taxon = levels[-1] + current_taxon
        taxons = [current_taxon]
        while current_taxon in taxon_to_parent:
            current_taxon = taxon_to_parent[current_taxon]
            taxons.append(current_taxon)
        taxons = list(reversed(taxons))

        print('\t'.join([id2, str('; '.join(taxons))]))
