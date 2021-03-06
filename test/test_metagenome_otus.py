#!/usr/bin/env python

import unittest
import subprocess
import os.path
import tempfile
import logging
import sys

import metagenome_otus

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','metagenome_otus.py')
logging.basicConfig(level=logging.DEBUG)

class Tests(unittest.TestCase):

#     def test_hello_world_script(self):
#         aln = '''>1_1_1_1
# KKK--
# >2_1_1_2
# --KFK
# >3_1_1_3
# --KK-'''
#         nucleotides = '''>1
# AAAAAAAAA
# >2
# AAATTTAAA
# >3
# AAAAAA'''
#         with tempfile.NamedTemporaryFile(suffix='aln.fa') as align_file:
#             align_file.write(aln)
#             align_file.flush()
#             with tempfile.NamedTemporaryFile(suffix='nuc.fa') as nuc_file:
#                 nuc_file.write(nucleotides)
#                 nuc_file.flush()
#                 
#                 cmd = "%s --alignment %s --reads %s --window_size 2" % (path_to_script,
#                                                                         align_file.name,
#                                                                         nuc_file.name)
#                 seqs = subprocess.check_output(cmd, shell=True)
#                 self.assertEqual(['AAATTT','AAAAAA'], seqs.split("\n"))
                
    def test_hello_world(self):
        aligned_sequences = []
        aligned_sequences.append(metagenome_otus.Sequence('1_1_1_1','KKK--'))
        aligned_sequences.append(metagenome_otus.Sequence('2_1_1_2','--KFK'))
        aligned_sequences.append(metagenome_otus.Sequence('3_1_1_3','--KK-'))
        nucs = {}
        nucs['1'] ='AAAAAAAAA'
        nucs['2'] ='AAATTTAAA'
        nucs['3'] ='AAAAAA'
        seqs = metagenome_otus.MetagenomeOtuFinder().find_windowed_sequences(aligned_sequences,
                                                                      nucs,
                                                                      2
                                                                      )
        self.assertEqual(['AAATTT','AAAAAA'], seqs)
        
    def test_gap_in_middle(self):
        aligned_sequences = []
        aligned_sequences.append(metagenome_otus.Sequence('1_1_1_1','KKK--'))
        aligned_sequences.append(metagenome_otus.Sequence('2_1_1_2','--KFK'))
        aligned_sequences.append(metagenome_otus.Sequence('3_1_1_3','--K-K'))
        nucs = {}
        nucs['1'] ='AAAAAAAAA'
        nucs['2'] ='AAATTTAAA'
        nucs['3'] ='AAAAAA'
        seqs = metagenome_otus.MetagenomeOtuFinder().find_windowed_sequences(aligned_sequences,
                                                                      nucs,
                                                                      3
                                                                      )
        self.assertEqual(['AAATTTAAA','AAA---AAA'], seqs)
        
    def test_specified_best_position(self):
        aligned_sequences = []
        aligned_sequences.append(metagenome_otus.Sequence('1_1_1_1','KKK--'))
        aligned_sequences.append(metagenome_otus.Sequence('2_1_1_2','--KFK'))
        aligned_sequences.append(metagenome_otus.Sequence('3_1_1_3','--K-K'))
        nucs = {}
        nucs['1'] ='AAAAAAAAA'
        nucs['2'] ='AAATTTAAA'
        nucs['3'] ='AAAAAA'
        seqs = metagenome_otus.MetagenomeOtuFinder().find_windowed_sequences(aligned_sequences,
                                                                      nucs,
                                                                      2,
                                                                      1
                                                                      )
        self.assertEqual(['AAAAAA'], seqs)
        
    def test_lower_case_insert(self):
        aligned_sequences = []
        aligned_sequences.append(metagenome_otus.Sequence('1_1_1_1','KKK---'))
        aligned_sequences.append(metagenome_otus.Sequence('2_1_1_2','--KFhK'))
        aligned_sequences.append(metagenome_otus.Sequence('3_1_1_3','--K--K'))
        nucs = {}
        nucs['1'] ='AAAAAAAAA'
        nucs['2'] ='AAATTTCACAAA'
        nucs['3'] ='AAAAAA'
        seqs = metagenome_otus.MetagenomeOtuFinder().find_windowed_sequences(aligned_sequences,
                                                                      nucs,
                                                                      3
                                                                      )
        self.assertEqual(['AAATTTAAA','AAA---AAA'], seqs)
        
    def test_lower_case_before_start_insert(self):
        aligned_sequences = []
        aligned_sequences.append(metagenome_otus.Sequence('1_1_1_1','KhKK---'))
        aligned_sequences.append(metagenome_otus.Sequence('2_1_1_2','-h-KFhK'))
        aligned_sequences.append(metagenome_otus.Sequence('3_1_1_3','---K--K'))
        nucs = {}
        nucs['1'] ='AAACACAAAAAA'
        nucs['2'] ='CACAAATTTCACAAA'
        nucs['3'] ='AAAAAA'
        seqs = metagenome_otus.MetagenomeOtuFinder().find_windowed_sequences(aligned_sequences,
                                                                      nucs,
                                                                      3
                                                                      )
        self.assertEqual(['AAATTTAAA','AAA---AAA'], seqs)
        seqs = metagenome_otus.MetagenomeOtuFinder().find_windowed_sequences(aligned_sequences,
                                                                      nucs,
                                                                      3,
                                                                      2
                                                                      )
        self.assertEqual(['AAATTTAAA','AAA---AAA'], seqs)
        
class TestSequence(unittest.TestCase):
    def test_orfm_nucleotides_revcom(self):
        s = metagenome_otus.Sequence('seq_2_5_1', 'LE')
        self.assertEqual('GGCATC', s.orfm_nucleotides('CGATGCC'))
        
if __name__ == "__main__":
    unittest.main()
