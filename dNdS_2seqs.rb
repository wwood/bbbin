#!/usr/bin/env ruby

# Given a fasta file with 2 dna sequences in it, work out the number of synonymous to non-synonymous mutations. Very simple!



require 'rubygems'
require 'bio'

seqs = []
    Bio::FlatFile.open(ARGF).each do |s|
      seqs.push s.seq
end

raise unless seqs.length ==2
raise unless seqs.collect{|s| s.length}.uniq.length==1

aln = Bio::Alignment.new(seqs)
	dN = 0
	dS = 0
      aln.each_window(3,3) do |seq_array|

	translates = seq_array.collect{|s| Bio::Sequence::NA.new(s).translate}
	if translates[0] == translates[1] and seq_array[0] != seq_array[1]
	dS += 1
elsif translates[0] != translates[1]
	dN += 1
end	
end  

puts ['dN','dS','dN/dS'].join("\t")
puts [dN,dS,dN.to_f/dS.to_f].join("\t")


