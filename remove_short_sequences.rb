#!/usr/bin/env ruby

# Sequences that came out of assembly sometimes are given the sequence
# SEQUENCETOOSHORTTOBLASTBPAFTERQUALITYCLIPPING
# This script removes those sequences from a fasta file

require 'bio'

Bio::FlatFile.open(ARGF).each do |seq|
  unless seq.seq.match(/SEQUENCETOOSHORTTOBLAST/i)
    puts seq.to_s
  end
end
