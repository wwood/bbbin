#!/usr/bin/env ruby

# Given a tblastn report from NCBI blast downloaded as a text Hit Table(text),
# and the fasta file of the nucleotide sequences, translate each of the
# sequences in the correct frame and print out that fasta file

require 'bio'
require 'faster_csv'

unless ARGV.length == 2
  $stderr.puts "Usage: #{$0} <hit table(text)> <fasta>"
  $stderr.puts "both files should be in the same order."
  exit
end

fastas = Bio::FlatFile.open(ARGV[1]).entries

i = 0
FasterCSV.foreach(ARGV[0], :col_sep => "\t") do |row|
  #skip comments
  next if row[0].nil? or row[0].match(/^\#/)

  # read the sequence name
  name = row[1]

  # read the sequence frame
  frame = row[4].split('/')[1].to_i

  # translate the sequence and print
  puts ">#{name}"
  puts Bio::Sequence::NA.new(fastas[i].seq).translate(frame)

  # finish iteration
  i += 1
  break if i>fastas.length
end
