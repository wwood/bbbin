#!/usr/bin/env ruby

require 'bio'

# Takes an nucleic acid sequence and chops it into codons

RARE_CODONS = {
  'AGA' => 'Arginine',
  'AGG' => 'Arginine',
  'ATA' => 'Isoleucine',
  'GGA' => 'Glycine',
  'CCC' => 'Proline',
  'CTA' => 'Leucine'
}

found = {}
all_seq = ARGF.read

Bio::Sequence::AA.new(all_seq).window_search(3, 3) do |window|
  seq = window.seq
  if RARE_CODONS.key?(window.seq)
    found[seq] ||= 0
    found[seq] += 1
  end
end

if found.length > 0
  found.each do |key, value|
    puts [key, RARE_CODONS[key], value].join(' ')
  end
else
  puts 'No rare codons found!'
end

puts
puts "Splits sequence:"
Bio::Sequence::AA.new(all_seq).window_search(3, 3) do |window|
  print "#{window.seq} "
end
puts
