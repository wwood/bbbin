#!/usr/bin/env python3

###############################################################################
#
#    Copyright (C) 2021 Ben Woodcroft
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
__copyright__ = "Copyright 2021"
__credits__ = ["Ben Woodcroft"]
__license__ = "GPL3"
__maintainer__ = "Ben Woodcroft"
__email__ = "benjwoodcroft near gmail.com"
__status__ = "Development"
__version__ = '0.0.0-dev'

import argparse
import logging
import sys
import os
import re

import pandas
from ete3 import Tree

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('-m','--metadata',help='GTDB metadata file',required=True)
    parent_parser.add_argument('-t','--tree', help='tree where tips are renamed', required=True)

    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    logging.info("{} version {}".format(os.path.basename(__file__), __version__))

    # Read metadata file
    df = pandas.read_csv(args.metadata, sep='\t')
    
    # Read tree file
    tree_str = open(args.tree).read()

    # Replace [] as ete cannot handle them
    # r = re.compile('\[.*?\]')
    # tree_str2 = r.sub('', tree_str)

    tree = Tree(tree_str)
    
    remove_start = re.compile('^.._')
    remove_end = re.compile('\.\d$')
    leaves = tree.get_leaves()
    for (accession, tax) in zip(df.accession, df.gtdb_taxonomy):
        a2 = remove_end.sub('',remove_start.sub('',accession))
        for l in leaves:
            if a2 in l.name:
                replacement = tax.replace(';','_')
                logging.debug("Replacing {} with {}, as it matched {}".format(l.name, replacement, a2))
                l.name = replacement

    print(tree.write())
    # import IPython; IPython.embed()
