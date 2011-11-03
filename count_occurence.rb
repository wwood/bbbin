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
# a - aligned - report occurrences within the alignment, column-wise
options = ARGV.getopts('d','i','a')

if options.values.select{|truthiness| truthiness}.length > 1
  $stderr.puts "Cannot use more than one of -i, -d and -a options together (at the moment), sorry"
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

alignment_hash = [] # for -a. array of hashes, where each array element is a position, and each hash is the number of occurrences of each AA at that position

Bio::FlatFile.foreach(ARGF) do |entry|
  last = nil
  seq = entry.seq
  i = 0
  seq.each_char do |char|
    if options['d'] # d for dipeptides
      if last
        i = "#{last}#{char}"
        hash[i] ||= 0
        hash[i] += 1
      end
      last = char
    elsif options['a']
      alignment_hash[i] ||= {}
      alignment_hash[i][char] ||= 0
      alignment_hash[i][char] += 1
    else
      hash[char] ||= 0
      hash[char] += 1
    end
    i += 1
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

if options['i']
  # do nothing
elsif options['a']
  puts amino_acids.join(',')
  alignment_hash.each do |hash2|
    amino_acids.each_with_index do |aa, i|
      print ',' unless i==0 #print comma before each count, except if this is the first count
      if hash2[aa]
        print hash2[aa]
      else
        print 0
      end
    end
    puts
  end
elsif options['d']
  pp hash
else
  amino_acids.each_with_index do |aa, i|
    print ',' unless i==0 #print comma before each count, except if this is the first count
    if hash[aa]
      print hash[aa]
    else
      print 0
    end
  end
  puts
  
  unfounds = hash.keys.reject{|a| amino_acids.include?(a)}
  if unfounds.length > 0
    $stderr.puts "Not reporting these amino acid occurrences: #{unfounds.join(', ')}."
  end
end