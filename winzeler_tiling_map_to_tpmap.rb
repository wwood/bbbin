#!/usr/bin/env ruby

# Convert the Winzeler Pftiling_tile-3.bpmap.txt file which is not a
# bpmap file into a tpmap file, which can then be converted into a
# bpmap file using the affy tiling array utils

# Taken from http://www.affymetrix.com/support/developer/powertools/changelog/gcos-agcc/tpmap.html

require 'rubygems'
require 'fastercsv'

#seq_group_name group_name
#version version
puts "#seq_group_name PfalciparumTiling"
puts "#version ben1"

index = 0
ARGF.each_line do |line|
  # skip header
  if index == 0
    index += 1
    next
  end
  row = line.strip.split("\t")
  
  raise unless row.length == 7

  seq = row[6]
  x = row[0]
  y = row[1]

  # sequences all have to be 25 apparently
  next unless seq.length == 25 # ignore coz I don't think the tpmap2bpmap (and beyond?) can handle them
  #  while seq.length < 25
  #    seq += "X"
  #  end
  raise unless seq.length == 25

  puts [
    seq,
    'f',
    "seq#{index}",
    index,
    x,
    y
  ].join(' ')
  
  index += 1
end
