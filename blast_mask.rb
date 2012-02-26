#!/usr/bin/env ruby

require 'rubygems'
require 'bio'
require 'csv'

# A script to mask out certain sections of a fasta file as
# specified by their homology to some blast database

class BlastHitArray < Array
  # Taking a string sequence, apply the mask encoded in the Array superclass
  # and return the masked string sequence
  def masked_sequence(sequence, surrounds = nil)
    seq = ''
    
     (1..sequence.length).each do |i|
      # run the gauntlet of masks, innocent until proven guilty
      masked = false
      each do |hit|
        if hit.contains?(i, surrounds)
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
  
  def contains?(one_based_index, surrounds = nil)
    surrounds ||= 0
    # If hit is a forward hit
    if from < to
      one_based_index >= from-surrounds and one_based_index <= to+surrounds
    else #If hit is a reverse hit
      one_based_index >= to-surrounds and one_based_index <= from+surrounds
    end
  end
end

if __FILE__ == $0
  require 'optparse'
  
  # Parse cmd line options
  USAGE = "Usage: blast_mask.rb [-m] <fasta_filename> <blast_filename>"
  options = {
  :print_masked_sequences_only => false,
  }
  o = OptionParser.new do |opts|
    opts.banner = USAGE
    
    opts.on("-m", "--masked-only", "Print out only those sequences that have blast hits (or equivalently, only those that are masked") do |v|
      options[:print_masked_sequences_only] = true
    end
    
    opts.on("-s", "--mask-surrounds AMOUNT_SURROUNDING", "Print out only those sequences that have blast hits (or equivalently, only those that are masked") do |v|
      options[:surrounding_to_mask] = v.to_i
      raise if options[:surrounding_to_mask] < 0
    end
  end
  o.parse!
  
  unless ARGV.length == 2
    $stderr.puts USAGE
    exit 1
  end
  
  fasta_filename = ARGV[0]
  blast_filename = ARGV[1]
  
  
  
  # read in blast data
  blasts = {} #hash of sequence identifiers to arrays of blasthits
  CSV.foreach(blast_filename, :col_sep => "\t") do |row|
    query = row[0]
    from = row[6].to_i
    to = row[7].to_i
    
    blasts[query] ||= BlastHitArray.new
    blasts[query].push Hit.new(from, to)
  end
  
  # read and print fasta file with the masking done
  Bio::FlatFile.foreach(Bio::FastaFormat, fasta_filename) do |entry|
    # Apply masking if required
    if blasts[entry.definition]
      masked = blasts[entry.definition].masked_sequence(entry.seq, options[:surrounding_to_mask])
      puts ">#{entry.definition}"
      puts masked
    else
      # No matching blast hits recorded. Print unless
      # the config says not to print these ones
      puts entry unless options[:print_masked_sequences_only]
    end
  end
end
