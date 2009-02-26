#!/usr/bin/env ruby

# Convert the Winzeler Pftiling_tile-3.bpmap.txt file which is not a
# bpmap file into a tpmap file, which can then be converted into a
# bpmap file using the affy tiling array utils

# Taken from http://www.affymetrix.com/support/developer/powertools/changelog/gcos-agcc/tpmap.html

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
  next if seq.length == 0

  # bump the index if there is any sequence, because that is what
  # sequenceNamer.pl is doing, and we want to line up with that
  index += 1
  i = index - 1

  # ignore seqs with length != 25 coz I don't think the tpmap2bpmap (and beyond?) can handle them
  # Cannot just add X's either because tpmap2bpmap can't handle those (tested)
  next unless seq.length == 25

  puts [
    seq,
    'f',
    "seq#{i}",
    i,
    x,
    y
  ].join(' ')
end
