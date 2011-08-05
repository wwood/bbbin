#!/usr/bin/env ruby

# Searching for silaffin-like proteins
# Modelled after Scheffel et. al. 2010 PNAS, with slight modifications
#
# A protein is considered a hit if it has a 100-2000 amino acid region that
# contains >=18% serine and >=10% lysine.
#
# Modifications:
# * Signal peptide is cleaved before window calculations are done, so that amino acids in the signal peptide cannot contribute to a window that passes the filter
# * SignalP yes/no is NN _or_ HMM, not both
# * SignalP cutoffs for NN and HMM are default
#
# Takes in a fasta file

require 'rubygems'
require 'bio'
require 'bio-signalp'

module Bio
  class SilaffinScreener
    # A protein is considered a hit if it has a 100-2000 amino acid region that
    # contains >=18% serine and >=10% lysine.
    def hit?(amino_acid_sequence_object)
      raise Exception, "Expected Bio::Sequence::AA object, found #{amino_acid_sequence_object.class}" unless amino_acid_sequence_object.kind_of?(Bio::Sequence::AA)
      
      # First, make sure that a signal peptide has been predicted
      result = Bio::SignalP::Wrapper.new.calculate(amino_acid_sequence_object.seq)
      return false unless result.signal?
      aa = Bio::Sequence::AA.new result.cleave(amino_acid_sequence_object.seq)
      
      # Second, make sure that there is a window size that fits the bill
      (100..2000).each do |window_size|
        aa.window_search(window_size) do |aa|
          if window_is_hit?(aa)
            return true
          end
        end
      end
      return false
    end

    # Takes a region of interest and see if that whole region is >=18% serine and
    # >= 10% lysine. Does NOT do the window search - use hit? for this (ie use
    # hit? for screening a whole protein)
    def window_is_hit?(aasequence)
      composition = aasequence.composition
      if composition['S'].to_f/aasequence.length >= 0.18 and composition['K'].to_f/aasequence.length >= 0.1
      return true
      else
      return false
      end
    end
  end
end

if $0 == __FILE__
  screener = Bio::SilaffinScreener.new
  
  unless ARGV.length == 1
    USAGE = 'Usage: silica_sliding_window.rb <fasta_filename>'
    $stderr.puts USAGE
    exit 1
  end

  Bio::FastaFormat.open(ARGV[0]).each do |sequence|
    hit = screener.hit?(sequence.aaseq)
    puts [
      sequence.entry_id,
      hit
    ].join("\t")
  end
end