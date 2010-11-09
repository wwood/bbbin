#!/usr/bin/env ruby

# Make each node in a tree have a unique name - otherwise FigTree at least
# cannot handle it

require 'bio'
require File.dirname(__FILE__) + '/uniq'

entries = Bio::FlatFile.open(ARGF).entries
uniq = Uniq.new

entries.each do |entry|
  entry.definition = uniq.make_unique(entry.definition);
  puts entry.to_s
end
