#!/usr/bin/env ruby

require 'rubygems'
require 'bio'
require 'fastercsv'

# A script to mask out certain sections of a fasta file as
# specified by their homology to some blast database

class BlastHitArray < Array
  # Taking a string sequence, apply the mask encoded in the Array superclass
  # and return the masked string sequence
  def masked_sequence(sequence)
    seq = ''
    
     (1..sequence.length).each do |i|
      # run the gauntlet of masks, innocent until proven guilty
      masked = false
      each do |hit|
        if hit.contains?(i)
          masked = true
          break
        end
      end
      
      # print a mask character or not, depending on verdict
      if masked
        seq += 'X'
      else
        seq += sequence[i-1..i-1] #sequence is 0-based index, whereas hits are 1-based
      end
    end
    return seq
  end
end

class Hit
  attr_accessor :from, :to
  def initialize(from, to)
    raise Exception, "Integers needed" unless from.kind_of?(Integer) and to.kind_of?(Integer)
    @from = from
    @to = to
  end
  
  def contains?(one_based_index)
    one_based_index >= from and one_based_index <= to
  end
end

if __FILE__ == $0
  fasta_filename = ARGV[0]
  blast_filename = ARGV[1]
  
  # read in blast data
  blasts = {} #hash of sequence identifiers to arrays of blasthits
  FasterCSV.foreach(blast_filename, :col_sep => "\t") do |row|
    query = row[0]
    from = row[6].to_i
    to = row[7].to_i
    
    blasts[query] ||= BlastHitArray.new
    blasts[query].push Hit.new(from, to)
  end
  
  # read and print fasta file with the masking done
  Bio::FlatFile.foreach(Bio::FastaFormat, fasta_filename) do |entry|
    # Apply masking if required
    if blasts[entry.entry_id]
      masked = blasts[entry.entry_id].masked_sequence(entry.seq)
      puts ">#{entry.entry_id}"
      puts masked
    else
      # No matching blast hits recorded
      puts entry
    end
  end
end