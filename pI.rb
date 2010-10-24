#!/usr/bin/env ruby

# Takes an input fasta file, and outputs the pI for each input protein

require 'rubygems'
require 'bio'
require 'isoelectric_point' #currently been refactored to fit as a bioruby plugin. 



Bio::FlatFile.foreach(Bio::FastaFormat, ARGF) do |entry|
  puts [
  entry.entry_id,
  entry.aaseq.calculate_iep,
  ].join("\t")
end
