#!/usr/bin/env ruby

# Ensure that the names of the sequences in a fasta file do not
# contain unwanted characters, such as /

require 'bio'
Bio::FlatFile.open(ARGF).each do |seq|
  name = seq.definition.gsub(/[\|\/\\]/, '_')
  puts ">#{name}"
  puts seq.seq
end
