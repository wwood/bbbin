#!/usr/bin/env ruby
# 
# A script to take a cazy list of proteins, and return a fasta file
# representing them.

require 'rubygems'
require 'fastercsv'
require 'bio'

serv = Bio::Fetch.new

FasterCSV.foreach("/home/ben/phd/cbm48/eukaryotes.cazy.csv", :col_sep => "\t", :headers => false) do |row|
  # take the first uniprot sequence, or failing that, the first genbank one. Maybe not a
  # perfect scenario, but eh for now
  
  seq = nil
  unless row[4].nil? or row[4].length == 2
    uniprot_id =  row[4].strip.split("\n")[0]
    
    system("wget http://www.uniprot.org/uniprot/#{uniprot_id}.fasta") #ain't REST beautiful
  else
    genbank_id = row[3].strip.split("\n")[0]
#    puts genbank_id
  end
end