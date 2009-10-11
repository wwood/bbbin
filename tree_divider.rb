#!/usr/bin/env ruby

# Input a tree, and divide each of the bootstraps on it by some number
# (10 by default), and then output the tree with a rounded integer number of
# that. So the point is to convert bootstrap numbers into percentages.

options = {
  :divisor => 10
}

# "Usage: tree_divider.rb <newick_tree_filename>"

require 'bio'

tree = Bio::FlatFile.open(Bio::Newick, ARGF).entries[0].tree

tree.each_edge do |node1, node2, edge|
  edge.distance = (edge.distance/options[:divisor]).round
end

puts tree.output(:newick)