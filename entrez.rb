#!/usr/bin/python
#
# Script to read FASTA sequences from GenBank
#
# Usage: entrez.rb   database  search_term output_file
#
# eg  entrez.rb protein nematode nematode.fna
#
# Note: user needs to place email address in the script - name@domain.edu.au, line 30
#
# Ben J Woodcroft 2009. Adapted from entrez.py by Ross Hall
#

require 'bio'
require 'reach'

puts Bio::NCBI::REST::EFetch.est(File.open(ARGV[0]).to_a.reach.strip!.to_i.retract, "fasta")

