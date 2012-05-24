#!/usr/bin/env ruby

# Take a facil file <fasta>.codons.txt which is a table of codon counts to consensus counts, and tally each of them. Then print the tally as a matrix

require 'csv'
require 'bio'

counts = {}
CSV.foreach(ARGV[0], :col_sep => "\t", :headers => true) do |row|
  codon = row[0]
  aa = row[8]
  counts[codon] ||= {}
  counts[codon][aa] ||= 0
  counts[codon][aa] += 1
end

aas = counts.values.collect {|aa_count| aa_count.keys}.flatten.uniq.sort

# Print aa header
print "\t"
puts [
  aas,
  'max',
  'expected',
  'different?'
  ].flatten.join("\t")

counts.each do |codon, aa_counts|
  print codon
  aas.each do |aa|
    if aa_counts.key?(aa)
      print "\t#{aa_counts[aa]}"
    else
      print "\t0"
    end
  end
  
  print "\t"
  max_aa = aa_counts.max{|ac1,ac2| ac1[1]<=>ac2[1]}[0]
  print max_aa
  
  print "\t"
  expected = Bio::Sequence::NA.new(codon).translate
  print expected
  
  print "\t"
  if max_aa=expected
    print "different"
  end
  puts
end
