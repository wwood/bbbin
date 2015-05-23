#!/usr/bin/env python

__author__ = "Ben Woodcroft"
__copyright__ = "Copyright 2015"
__credits__ = ["Ben Woodcroft"]
__license__ = "GPL3"
__maintainer__ = "Ben Woodcroft"
__email__ = "b.woodcroft near uq.edu.au"
__status__ = "Development"

import argparse
import sys
import os
from Bio.Seq import Seq
import itertools
import logging
import re
import IPython

# Stolen from https://github.com/lh3/readfq/blob/master/readfq.py
def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break
            
            
def find_best_window(protein_alignment, stretch_length):
    # Convert the alignment into a True/False matrix for ease,
    # True meaning that there is something aligned, else False
    binary_alignment = []
    for s in protein_alignment:
        aln = []
        for i, base in enumerate(s.seq):
            if base=='-':
                aln.append(False)
            else:
                aln.append(True)
        binary_alignment.append(aln)
    print binary_alignment
    
    # Find the number of aligned bases at each position
    current_best_position = 0
    current_max_num_aligned_bases = 0
    for i in range(0, len(binary_alignment[0])-stretch_length+1):
        num_bases_covered_here = 0
        end_index = i+stretch_length-1 
        for s in binary_alignment:
            if not s[i] or not s[end_index]: continue #ignore reads that don't cover the entirety
            num_bases_covered_here += sum(s[i:end_index])
        logging.debug("Found %i aligned bases at position %i" % (num_bases_covered_here, i))
        if num_bases_covered_here > current_max_num_aligned_bases:
            current_best_position = i
            current_max_num_aligned_bases = num_bases_covered_here
    logging.info("Found a window starting at position %i with %i bases aligned" % (current_best_position,
                                                                                   current_max_num_aligned_bases
                                                                                   ))
    return current_best_position
            
            
class Sequence:
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq
        
    def unaligned_length(self):
        return len(re.sub('-','',self.seq))
        
    def un_orfm_name(self):
        return re.sub('_\d+_\d+_\d+$', '', self.name)
    
    def orfm_nucleotides(self, original_nucleotide_sequence):
        m = re.search('_(\d+)_(\d+)_\d+$', self.name)
        start = int(m.groups(0)[0])-1
        translated_seq = original_nucleotide_sequence[start:(start+3*self.unaligned_length())]
        logging.debug("Returning orfm nucleotides %s" % translated_seq)
        if int(m.groups(0)[1]) > 3:
            # revcomp type frame
            return(str(Seq(self.seq).reverse_complement()))
        else:
            return(translated_seq)
        
def nucleotide_alignment(protein_sequence, nucleotides, start_position, stretch_length):
    '''Line up the nucleotides and the proteins, and return the alignment 
    across start_position for stretch length amino acids'''
    
    codons = []
    # For each position in the amino acid sequence
    # If non-dash character, take 3 nucleotides off the nucleotide sequence and
    # add that as the codon
    # else add None
    for aa in protein_sequence.seq:
        if aa=='-':
            codons.append(None)
        else:
            if len(nucleotides) < 3: raise Exception("Insufficient nucleotide length found")
            codons.append(nucleotides[:2])
            nucleotides = nucleotides[:2]
    if len(nucleotides) > 0: raise Exception("Insufficient protein length found")
    
    return ''.join(itertools.chain(codons[start_position:(start_position+stretch_length-1)]))
    
class MetagenomeOtuFinder:
    def __init__(self):
        pass
    
    def find_windowed_sequences(self,
                                aligned_sequences,
                                nucleotide_sequences,
                                stretch_length
                                ):
        # Find the stretch of the protein that has the most number of aligned bases in a 20 position stretch,
        # excluding sequences that do not aligned to the first and last bases
        best_position = find_best_window(aligned_sequences, stretch_length)
        logging.info("Found best section of the alignment starting from %i" % (best_position+1))
        
        # For each read aligned to that region i.e. has the first and last bases,
        # print out the corresponding nucleotide sequence
        windowed_sequences = []
        for s in aligned_sequences:
            if s.seq[best_position] != '-' and s.seq[best_position+stretch_length-1] != '-':
                nuc = nucleotide_sequences[s.un_orfm_name()]
                windowed_sequences.append(
                                          nucleotide_alignment(s, s.orfm_nucleotides(nuc), best_position, stretch_length)
                                          )
        return windowed_sequences

if __name__ == '__main__':
    parser = argparse.ArgumentParser(add_help=False)
    
    parser.add_argument('--alignment', metavar='aligned_fasta', help="Protein sequences hmmaligned and converted to fasta format with seqmagick", required=True)
    parser.add_argument('--reads', metavar='raw_reads', help='Unaligned nucleotide sequences that were translated into the protein sequences', required=True)
    parser.add_argument('--window_size', metavar='bp', help='Number of base pairs to use in continuous window', default=20)
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG)
    
    # Read in the fasta Alignment
    protein_alignment = []
    for name, seq, _ in readfq(args.alignment):
        protein_alignment.append(Sequence(name, seq))
    logging.info("Read in %i aligned protein sequences e.g. %s %s" % (len(protein_alignment),
                                                          protein_alignment[0].name,
                                                          protein_alignment[0].seq
                                                          ))
    
    # Read in the original nucleotide sequences
    nucleotide_sequences = {}
    for name, seq, _ in readfq(args.reads):
        eg_name = name
        nucleotide_sequences[name] = seq
    logging.info("Read in %i nucleotide sequences e.g. %s %s" % (len(nucleotide_sequences),
                                                          eg_name,
                                                          nucleotide_sequences[eg_name]
                                                          ))
    
    aligned_sequences = find_windowed_sequences(protein_alignment,
                                                nucleotide_sequences,
                                                args.window_size)
    print '\n'.join(aligned_sequences)
    

    
    
    
    
    
    
    
    
    