#!/usr/bin/env python2.7

import dendropy

import argparse
import logging
import os
import sys

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', help='tree file to analyse', required=True)
    parser.add_argument('--clades', help='clades to find monophyly between', nargs='+', required=True)
    parser.add_argument('--separator', help="character separating genomes from gene in IDs eg '~'", required=True)
    args = parser.parse_args()

    separator = args.separator

    # Read each file of clades
    genome_to_clade = {}
    for clade in args.clades:
        with open(clade) as f:
            for genome in f.readlines():
                g = genome.replace('_',' ').strip()
                if g != '':
                    genome_to_clade[g+separator] = clade

    # Open tree
    tree = dendropy.Tree.get(
        path=args.t,
        schema='newick')

    def clade_of_node(node, genome_to_clade):
        l = node.taxon.label
        for genome in genome_to_clade.keys():
            g = node.taxon.label
            if g.startswith(genome):
                return genome_to_clade[genome]
        return None

    def monophyletic_one_clade(node, clade, genome_to_clade):
        for n in node.leaf_nodes():
            if n.taxon is not None:
                found_clade = clade_of_node(node, genome_to_clade)
                if found_clade is not None and found_clade != clade:
                    return False
        return True

    def clades_of_children(node, genome_to_clade):
        founds = set()
        for n in node.leaf_nodes():
            found_clade = clade_of_node(n, genome_to_clade)
            founds.add(found_clade)
        return founds

    # Find a node that is not part of any genome
    node_to_reroot_on = None
    for tip in tree.leaf_node_iter():
        if clade_of_node(tip, genome_to_clade) is None:
            node_to_reroot_on = tip
            break
    if node_to_reroot_on is None: # Entire tree is monophyletic, so stop here
        clades = set()
        for tip in tree.leaf_node_iter():
            clades.add(clade_of_node(tip, genome_to_clade))
        print '%s\t%s\t%s' % (args.t, ','.join(sorted(clades)), 'all')

    else:
        # Reroot on it
        tree.reroot_at_node(node_to_reroot_on)
        if node_to_reroot_on.parent_node != None: raise Exception

        printed_clades = set()

        # for each label
        for tip in tree.leaf_node_iter():
            # if the label is in one of the clades of interest
            genome = None
            label = tip.taxon.label
            for g in genome_to_clade.keys():
                if label.startswith(g):
                    genome = g

            if genome is not None:
                mrca = tip
                while mrca is not None and None not in clades_of_children(mrca, genome_to_clade):
                    last_mrca = mrca
                    mrca = mrca.parent_node

                final_clades = clades_of_children(last_mrca, genome_to_clade)
                leaf_names = ', '.join(sorted([n.taxon.label for n in last_mrca.leaf_nodes()]))
                if leaf_names not in printed_clades:
                    printed_clades.add(leaf_names)
                    if len(final_clades) > 1:
                        print "%s\t%s\t%s" % (
                            args.t,
                            ",".join(sorted(final_clades)),
                            leaf_names)
