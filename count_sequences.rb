#!/usr/bin/ruby

# Count the number of sequences in a fasta file

require 'bio'

seqs = Bio::FlatFile.open(nil, $stdin)
count = 0

for s in seqs
  count = count +1
end

puts count