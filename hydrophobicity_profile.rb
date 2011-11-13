#!/usr/bin/env ruby

require 'rubygems'
require 'bio'
require 'hydrophobicity'
require 'optparse'

options = {
  :hydrophibic_binary => false
}
OptionParser.new do |opts|
  opts.banner = "Usage: hydrophobicity_profile.rb [--hydrophibic-binary]\n\tFASTA file is either piped in or specified as an argument"

  opts.on(nil, "--hydrophibic-binary", "Do a simple yes/no is it hydrophobic? at each position") do |v|
    options[:hydrophibic_binary] = true
  end
end.parse!

# Print out a csv of the hydrophobicity of a protein sequence, fed in by fasta file.
h = Hydrophobicity.new
Bio::FlatFile.foreach(ARGF) do |entry|
  if options[:hydrophibic_binary]
    array = []
    entry.seq.each_char do |char|
      if %w(A I L M F W Y V).include?(char)
        array.push 1
      else
        array.push 0
      end
    end
    puts array.join(",")
  else
    puts h.hydrophobicity_profile(entry.seq).join(',')
  end
end
