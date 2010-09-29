#!/usr/bin/env ruby

# Count the total occurance of each letter (i.e. amino acid or nucleotide) in the
# input fasta file.

require 'rubygems'
require 'bio'
require 'pp'
require 'optparse'

hash = {}
# d - use dipeptides
# i - report occurences for each sequence individually
options = ARGV.getopts("d",'i')

if options['i'] and options['d']
  $stderr.puts "Cannot use -i and -d options together at the moment."
  exit 1
end

amino_acids = Bio::AminoAcid::Data::WEIGHT.keys
# Print headers for 1 sequence at a time
if options['i']
  puts [
  'Identifier',
  amino_acids
  ].flatten.join("\t")
end

Bio::FlatFile.foreach(ARGF) do |entry|
  last = nil
  seq = entry.seq
  seq.each_char do |char|
    if options['d'] # d for dipeptides
      if last
        i = "#{last}#{char}"
        hash[i] ||= 0
        hash[i] += 1
      end
      last = char
    else
      hash[char] ||= 0
      hash[char] += 1
    end
  end
  
  # if reporting for each seq, report now
  if options['i']
    if hash.empty?
      $stderr.puts "Ignoring #{entry.definition} since it appears to have zero length"
    else
      print entry.definition
      amino_acids.each do |aa|
        to_print = hash[aa]
        to_print = 0 if to_print.to_s=='' #print 0 not nothing to make later analyses easier 
        print "\t#{to_print}"
      end
      puts
    end
    hash = {} # reset the hash for the next sequence
  end
end

pp hash unless options['i']
