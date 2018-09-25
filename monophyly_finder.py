#!/usr/bin/env python2.7

import dendropy

import argparse
import logging
import os
import sys

debug={1: logging.CRITICAL,
       2: logging.ERROR,
       3: logging.WARNING,
       4: logging.INFO,
       5: logging.DEBUG}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', help='tree file to analyse', required=True)
    parser.add_argument('--clades', help='clades to find monophyly between', nargs='+', required=True)
    parser.add_argument('--separator', help="character separating genomes from gene in IDs eg '~'", required=True)
    parser.add_argument('--log', help='Output logging information to file', type=str, default=False)
    parser.add_argument('--verbosity', help='1 - 5, 1 being silent, 5 being noisy indeed. Default = 4', type=int, default=4)

    args = parser.parse_args()

    if args.log:
        if os.path.isfile(args.log):    
            raise Exception("File %s exists" % args.log)
        logging.basicConfig(filename=args.log, level=debug[args.verbosity], format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    else:
        logging.basicConfig(level=debug[args.verbosity], format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    
    logging.debug('Starting monophyly finder')

    separator = args.separator
    logging.debug('Separator is: %s' % args.separator)
    
    # Read each file of clades
    logging.debug('Reading in %i clades file(s)' % len(args.clades))
    genome_to_clade = {}
    for idx, clade in enumerate(args.clades):
        logging.debug('Reading clade file %i' % (idx+1))
        with open(clade) as f:
            for f_idx, genome in enumerate(f.readlines()):
                g = genome.replace('_',' ').strip()
                if g != '':
                    logging.debug('Found genome: %s' % (g))
                    genome_to_clade[g+separator] = clade
                else:
                    logging.debug('Found missing line in file at line %i' % f_idx)
    # Open tree
    logging.debug('Reading in tree')
    tree = dendropy.Tree.get(
        path=args.t,
        schema='newick')
    logging.debug('Read in tree with %i tips' % len(tree.leaf_nodes()))

    def clade_of_node(node, genome_to_clade):
        for genome in genome_to_clade.keys():
            g = node.taxon.label
            if g.startswith(genome):
                return genome_to_clade[genome]
        return None

    def clades_of_children(node, genome_to_clade):
        founds = set()
        for n in node.leaf_nodes():
            found_clade = clade_of_node(n, genome_to_clade)
            founds.add(found_clade)
        return founds

    # Find a node that is not part of any genome
    logging.debug('Finding if tree needs rerooting')
    node_to_reroot_on = None
    for tip in tree.leaf_node_iter():
        if clade_of_node(tip, genome_to_clade) is None:
            node_to_reroot_on = tip
            logging.debug('Tree needs rerooting on tip: %s' % tip.taxon.label)
            break

    if node_to_reroot_on is None: 
        logging.debug('Entire tree is monophyletic, no rerooting required')
        clades = set()
        for tip in tree.leaf_node_iter():
            clades.add(clade_of_node(tip, genome_to_clade))
        print '%s\t%s\t%s' % (args.t, ','.join(sorted(clades)), 'all')

    else:
        # Reroot on it
        logging.debug('Rerooting')
        tree.reroot_at_node(node_to_reroot_on)
        if node_to_reroot_on.parent_node != None: raise Exception

        printed_clades = set()

        # for each label
        logging.debug('Iterating through tree tips')
        for idx, tip in enumerate(tree.leaf_node_iter()):
            # if the label is in one of the clades of interest
            idx+=1
            logging.debug('On tip: %i' % idx)
            genome = None
            label = tip.taxon.label
            for g in genome_to_clade.keys():
                if label.startswith(g):
                    genome = g

            if genome is not None:
                logging.debug('FOUND: Tip in defined groups: %s' % genome)
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
            else:
                logging.debug('Tip not in defined groups')