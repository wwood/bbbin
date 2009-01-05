#!/usr/bin/env ruby
# 
# Take the data from the Slovakian Conference paper
# send to me by David Stapleton and return a fasta file of all of the CBMs

require 'rubygems'
require 'bio'
require 'fastercsv'

hash = Bio::FlatFile.open('/home/ben/phd/cbm48/7/SlovakianConferenceTableGenbanks.fa').to_hash(:accession)

FasterCSV.foreach('/home/ben/phd/cbm48/7/SlovakianConferenceTable.csv', :col_sep => ' ') do |row|
  row.reject!{|e| e.nil?}
  
  next if row.length != 8 #ignore category lines
  
  genbank = row[4]
  length = row[5].to_i
  start = row[6].to_i
  stop = row[7].to_i
  type = row[0]
  species = "#{row[2]} #{row[3]}"
  
  seq = hash[genbank]
  aas = seq.seq
  
  puts ">#{genbank} #{type} #{species}"
  puts aas[(start-1)..(stop-1)]
end
