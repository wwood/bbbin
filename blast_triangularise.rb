#!/usr/bin/env ruby

# Self vs. self blast results essentially give square 
# distance matrices, because each protein pair is 
# compared twice - once for each order. 
# E.g. protein1 vs. protein2
# and 
# protein2 vs. protein1.
#
# This script parses a blast -m 8 (tabular) result file
# and removes half of the square (making it a triangle),
# and additionally removes the diagonal (e.g. protein1 
# vs protein1) alignments.

require 'rubygems'
require 'fastercsv'

if __FILE__ == $0
  pairs = {}
  FasterCSV.foreach(ARGV[0], :col_sep => "\t") do |row|
    p1 = row[0]
    p2 = row[1]
    index = "#{p1} #{p2}"
    index_reverse = "#{p2} #{p1}"
    next if pairs[index_reverse] #remove 1 triangle
    next if p1 == p2 #remove diagonal
    pairs[index] = true
    puts row.join("\t")
  end
end