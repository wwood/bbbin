#!/usr/bin/env ruby

# trimpoly.rb - trimpoly wrapper
#
# Copyright::  Copyright (C) 2008  Ben J Woodcroft <donttrustben at gmail.com>
# 
# License::    The Ruby License

require 'tempfile'

module Bio
  class Trimpoly
    def self.trim_both(sequence)
      Tempfile.open('trimpoly1temp') do |file1|
        # write sequence to file
        file1.puts '>trimpolytemp'
        file1.puts sequence
        file1.flush
        Tempfile.open('trimpoly1temp') do |file2|
          result = system "trimpoly -Y <#{file1.path} >#{file2.path}"
          
          line = File.open(file2.path).readlines[0]
          
          splits = line.split("\t")
          start = splits[2].to_i
          stop = splits[3].to_i
          
          if stop > start
            return sequence[start-1..stop-1]
          else
            return nil
          end
        end
      end
    end
  end
end


if $0 == __FILE__
  require 'bio'
  
  Bio::FlatFile.open(ARGV[0]).each do |seq|
    trimmed = Bio::Trimpoly.trim_both(seq.seq)
    if trimmed
      puts ">#{seq.entry_id}"
      puts trimmed
    end
  end
end