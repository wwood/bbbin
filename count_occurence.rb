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
options = ARGV.getopts("d",'i:')
 
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
result = hash[options['i']]
if result.nil?
puts 0
else
puts result
end
hash = {}
end
 end
 
 pp hash unless options['i']
