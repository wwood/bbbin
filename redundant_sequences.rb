#!/usr/bin/env ruby
# 
# Find redundant sequences in a file, to weed out strange ones

require 'bio'
require 'rubygems'
require 'reach'
require 'optparse'

USAGE = [
  'Usage: redundant_sequences.rb [-p] <fasta>',
]
options = {
    :print_all => false,
}
o = OptionParser.new do |opts|
  opts.banner = USAGE
  
  opts.on('-p','--print-all', "print out the redundancy reduced sequences instead of the redundant ones, printing errors to STDOUT if the definitions aren't identical") do |v|
    options[:print_all] = true
  end
end
o.parse!

hash = {}

# Collect everything
Bio::FlatFile.open(ARGF).each do |seq|
  hash[seq.seq] ||= []
  hash[seq.seq].push seq
end

# Print out if there is duplicates
hash.each do |key, seqs|
  if options[:print_all]
    if seqs.reach.definition.uniq.length > 1
      $stderr.puts "Sequences identical but definitions not: `#{seqs.reach.definition.uniq.join('\', `')}'"
    end
    puts seqs[0]
  else
    if seqs.length > 1
      puts
      puts ">#{seqs.reach.definition.join("\n>")}"
      puts seqs[0].seq
    end
  end
end
