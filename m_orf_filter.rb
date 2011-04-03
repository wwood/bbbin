#!/usr/bin/env ruby

# A script to remove sequences from a FASTA that don't start with methionine

require 'bio'

if __FILE__ == $0
  Bio::FlatFile.foreach(ARGF) do |entry|
    if entry.seq.match(/^M/)
      puts entry
    end
  end
end