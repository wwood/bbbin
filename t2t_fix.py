#!/usr/bin/env python2.7

import logging
import os
import sys
import argparse
import re
from string import split as _
import IPython

parser = argparse.ArgumentParser(description='''--- t2t_fix %s --- fix taxonomy files so that they validate with tax2tree's validate mode''')
parser.add_argument('-i', '--input_taxonomy', help='output the taxonomy to this file', required=True)
parser.add_argument('-o', '--output_taxonomy', help='output the taxonomy to this file', required=True)
parser.add_argument('--debug', help='output debug information', action="store_true")

args = parser.parse_args()

if args.debug:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)


# read in taxonomy file, keeping hashes of taxonomy => [[parent,count],..] and just an array of the taxonomies
taxonomy = []
expected_taxonomy_levels = 7
rank_bits = _('k p c o f g s')

child_to_parent_counts = [None] * (expected_taxonomy_levels-1)

for i in range(expected_taxonomy_levels-1):
    child_to_parent_counts[i] = {}
    
regex = re.compile(r'^.__')
    
logging.info("Reading in taxonomy file..")
with open(args.input_taxonomy) as f:
    for input_line1 in f:
        input_line = input_line1.rstrip()
        # split first by tab, then by '; '
        splits = input_line.split("\t")
        if len(splits) != 2:
            logging.error("Found line that didn't have 2 tab-separated components: %s" % input_line)
            exit(1)
        splits2 = splits[1].split('; ')
        if len(splits2) != expected_taxonomy_levels:
            logging.error("Found line that didn't have the expected number of (%s) components: %s" % (expected_taxonomy_levels, input_line))
            exit(1)
            
        # get rid of extra bits at the beginning
        splits3 = []
        for s in splits2:
            reg = regex.match(s)
            if reg:
                splits3.append(s[reg.end():])
            elif len(s) > 0:
                logging.warn("Found unexpected form for taxonomy in %s", str(s))
                splits2.append(s)
            # if now empty, then don't get parents at all

            
        for i, lineage in enumerate(splits3):
            if i==0: continue #first levels don't have parents
            if lineage=='': continue #empty taxon strings are not relevant
            
            #stupid python, this would be two shorter lines in ruby: '||=' operator ftw
            parent = splits3[i-1]
            if not child_to_parent_counts[i-1].has_key(lineage):
                child_to_parent_counts[i-1][lineage] = {}
            if not child_to_parent_counts[i-1][lineage].has_key(parent):
                child_to_parent_counts[i-1][lineage][parent] = 0
            
            child_to_parent_counts[i-1][lineage][parent] += 1
            
        taxonomy.append([splits[0], splits3])

# work out all the parents by majority vote, except for the last entry
logging.info("Working out parents..")
parents = [None] * (expected_taxonomy_levels-2)
num_with_multiple_parents = 0
for i in range(expected_taxonomy_levels-2):
    parents[i] = {}
    for lineage, parents_hash in child_to_parent_counts[i].iteritems():
        max_parent = ''
        max_parent_count = 0
        if len(parents_hash) > 1:
            num_with_multiple_parents += 1
        for parent, count in parents_hash.iteritems():
            if count > max_parent_count:
                max_parent = parent
                max_parent_count = count
        parents[i][lineage] = max_parent
logging.info("Found %s taxons with multiple parents" % num_with_multiple_parents)
    
# go through each line, printing out the 'fixed' taxonomy
logging.info("Writing output file..")
with open(args.output_taxonomy,'w') as f:
    for entry in taxonomy:
        gg_id = entry[0]
        original_tax = entry[1]

        # iterate from highest resolution to least resolution
        to_print = [None]*expected_taxonomy_levels
        
        for rev_i, lineage in enumerate(reversed(original_tax)):
            i = expected_taxonomy_levels-rev_i-1
            if i==0:
                continue
            if i == expected_taxonomy_levels-1:
                # when printing the last column, change it so that genus and species both written, to resolve epithet conflicts
                genus = original_tax[i-1]
                if original_tax[i] == '':
                    to_print[i] = "s__"
                else:
                    to_print[i] = "s__%s %s" %(genus, original_tax[i])
                to_print[i-1] = "%s__%s" % (rank_bits[i-1], original_tax[i-1])
            else:
                # here's where conflicting names are searched for
                # TODO: warn if multiple parents change
                lineage = to_print[i][regex.match(to_print[i]).end():]
                if lineage == "":
                    to_print[i-1] = "%s__%s" % (rank_bits[i-1], original_tax[i-1])
                else:
                    parent = parents[i-1][lineage]
                    to_print[i-1] = "%s__%s" % (rank_bits[i-1], parent)
                    
        # Check if the kingdom, phylum or class levels have changed - these are 
        # possibly horribly damaging cases if not done perfectly.
        original_tax_with_qualifiers = [rank_bits[i] + '__' + tax for i, tax in enumerate(original_tax)]
        if to_print[:2] != original_tax_with_qualifiers[:2]:
            logging.warn("Watch out for this one - is the large taxonomic change from ")
            logging.warn(gg_id + "\t" + '; '.join(original_tax_with_qualifiers))
            logging.warn("to")
            logging.warn(gg_id + "\t" + '; '.join(to_print))
            logging.warn("sane?")
                
        f.write(gg_id)
        f.write("\t")
        f.write('; '.join(to_print))
        f.write("\n")
        
logging.info("Finished")
        
        
        