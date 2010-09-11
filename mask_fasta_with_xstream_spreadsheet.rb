#!/usr/bin/env ruby

# Takes a fasta file and a XSTREAM output spreadsheet.
# Prints out the fasta file with the repeats masked out with an X.

require 'rubygems'
require 'bio'
require 'fastercsv'

class Repeat
  attr_accessor :start, :stop
end

if __FILE__ == $0
  unless ARGV.length == 2
    $stderr.puts "Usage: mask_fasta_with_xstream_spreadsheet.rb <fasta_file> <xstream_output_spreadsheet>"
    exit
  end
  fasta_filename = ARGV[0]
  xstream_output_filename = ARGV[1]
  
  proteins = {} # hash of protein id to an array of Repeat objects
  
  # Read in each of the repeats from the xstream file
  FasterCSV.foreach(xstream_output_filename, :col_sep => "\t", :headers => true) do |row|
    id = row[0]
    start = row[2].to_i
    stop = row[3].to_i
    repeat = Repeat.new
    repeat.start = start
    repeat.stop = stop
    proteins[id] ||= []
    proteins[id].push repeat
  end
  
  # for each protein in the fasta file, output it masked
  Bio::FlatFile.open(Bio::FastaFormat, fasta_filename).each_entry do |f|
    if proteins[f.definition]
      # a repeat has been detected. Print with mask
      masked_sequence = f.seq
      proteins[f.definition].each do |repeat|
        (repeat.start..repeat.stop).each do |index|
          masked_sequence[index-1] = 'X'
        end
      end
      puts ">#{f.definition}"
      puts masked_sequence
    else
      puts f
    end
  end
end