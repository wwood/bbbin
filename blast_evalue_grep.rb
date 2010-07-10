#!/usr/bin/env ruby

# Takes in a blast -m 8 file, and then returns lines with a better
# than defined e-value

if __FILE__ == $0
  require 'rubygems'
  require 'optparse'
  
  # Parse cmd line options
  USAGE = "Usage: blast_mask.rb [-m] <fasta_filename> <blast_filename>"
  cutoff = 1e-5
  options = {
  :print_masked_sequences_only => false,
  }
  o = OptionParser.new do |opts|
    opts.banner = USAGE
    
    opts.on("-e", "--e-value [E-VALUE]", "Print out sequences with better (i.e. smaller) e-values in the file") do |v|
      cutoff = v.to_f
    end
  end
  o.parse!
  
  ARGF.each do |line|
    splits = line.split("\t")
    evalue = splits[10].to_f
    if evalue <= cutoff
      print line
    end
  end
end