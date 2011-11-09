#!/usr/bin/env ruby

# Take a nucleotide fasta file, and output the A/T content of each codon.

if __FILE__ == $0
  require 'rubygems'
  require 'bio'
  
  Bio::FlatFile.foreach(ARGF) do |entry|
    seq = entry.to_biosequence
    
    print entry.entry_id
    count = 0
    seq.window_search(3,3) do |codon|
      if codon.match(/^[ATGC]{3}$/)
        atgc = {}
        %w(A T G C).each do |char|
          atgc[char] = 0
        end
        codon.each_char do |c|
          atgc[c] += 1
        end
        print ",#{atgc['A']+atgc['T']}"
      else
        print ','
      end
      
      count += 1
      break if count > 25
    end
    
    puts
  end
end