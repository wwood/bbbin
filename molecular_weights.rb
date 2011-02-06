#!/usr/bin/env ruby

require 'rubygems'
require 'bio'
  
if __FILE__ == $0
  Bio::FlatFile.foreach(ARGF) do |seq|
    aa_seq = seq.seq.gsub(/\*.*/,'').gsub(/X/,'') #Otherwise there is Runtime Errors thrown sometimes and not all of the file is parsed
    puts [
    Bio::Sequence::AA.new(aa_seq).molecular_weight,
    seq.definition,
    ].join("\t")
  end
end