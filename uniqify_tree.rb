#!/usr/bin/env ruby

# Make each node in a tree have a unique name - otherwise FigTree at least
# cannot handle it

require 'bio'
require '/home/ben/bin/uniq'

tree = Bio::FlatFile.open(ARGF).entries[0].tree
uniq = Uniq.new

tree.each_node do |node|
  #I get internal nodes here - not sure how else to skip
  next if node.name.nil? or node.name.length == 0
  
  node.name = uniq.make_unique(node.name);
end
