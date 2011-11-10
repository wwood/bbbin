#!/usr/bin/env ruby

# Take a nucleotide fasta file, and output the A/T content of each codon.

if __FILE__ == $0
  require 'rubygems'
  require 'bio'
  require 'optparse'
  
  options = ARGV.getopts('n:hl:')
  number = options['n'].to_i
  number ||= 3 #default to codons
  
  length = options['l'].to_i
  if length.nil? or length == 0
    length = 25
  end
  
  if options['h']
    $stderr.puts "Usage: at_bias.rb [-n <number_of_bases_to_search>] [-l length] <fasta_path>"
    exit 1
  end
  
  Bio::FlatFile.foreach(ARGF) do |entry|
    seq = entry.to_biosequence
    
    print entry.entry_id
    count = 0
    seq.window_search(number,number) do |codon|
      atgc = {}
      %w(A T G C).each do |char|
        atgc[char] = 0
      end
      
      if codon.match(/^[ATGC]{#{number}}$/)
        codon.each_char do |c|
          atgc[c] += 1
        end
      end

      print ",#{atgc['A']+atgc['T']}"
      
      count += 1
      break if count > length
    end
    
    puts
  end
end