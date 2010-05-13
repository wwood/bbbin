#!/usr/bin/env ruby

# Count the total occurance of each letter (i.e. amino acid or nucleotide) in the
# input fasta file.

 require 'rubygems'
 require 'bio'
 require 'pp'
 
 hash = {}
 
 Bio::FlatFile.foreach(ARGF) do |entry|
   entry.seq.each_char do |char|
     hash[char] ||= 0
     hash[char] += 1
   end
 end
 
 pp hash