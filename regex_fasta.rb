#!/usr/bin/env ruby

# Pretty simple idea - given a regex and a fasta file, return the sequences of all those sequences that match the regex


require 'rubygems'
require 'bio'
require 'optparse'



if __FILE__ == $0
  options = ARGV.getopts('r:n')
  regex = options['r']
  unless regex
    $stderr.puts "Usage: regex_fasta.rb -r <regex> <fasta_path>"
    exit 1
  end
  
  Bio::FlatFile.foreach(ARGF) do |entry|
    if options['n']
      to_print = [entry.definition]
      if matches = entry.seq.match(regex)
        to_print.push $`.length+1
      else
        to_print.push '-1' 
      end
      puts to_print.join("\t")
    else
      if entry.seq.match(regex)
        puts entry
      end
    end
  end
end