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

    def self.trim_poly_a(sequence)
      Tempfile.open('trimpoly1temp') do |file1|
        # write sequence to file
        file1.puts '>trimpolytemp'
        file1.puts sequence
        file1.flush
        Tempfile.open('trimpoly1temp') do |file2|
          result = system "trimpoly -e 3 <#{file1.path} >#{file2.path}"

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

    def self.trim_poly_t(sequence)
      Tempfile.open('trimpoly1temp') do |file1|
        # write sequence to file
        file1.puts '>trimpolytemp'
        file1.puts sequence
        file1.flush
        Tempfile.open('trimpoly1temp') do |file2|
          result = system "trimpoly -e 5 <#{file1.path} >#{file2.path}"

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
  require 'optparse'

  trim_method = nil

  USAGE = 'Usage: trimpoly.rb [-AT] <fasta_file>'
  OptionParser.new do |opts|
    opts.banner = USAGE

    opts.on('-A', "--trim_polyA", "Only trim the 3' end of polyA. Do not trim the 5' end of polyT (default is to do both)") do
      trim_method = [:trim_poly_a]
    end

    opts.on('-T', "--trim_polyT", "Only trim the 5' end of polyT. Do not trim the 3' end of polyA (default is to do both)") do
      raise OptionParser::ParseError, "Only one of -T, -A and -Y are allowed." unless trim_method.nil?
      trim_method = [:trim_poly_t]
    end

    opts.on('-Y', "--trim_both", "Trim both polyA and polyT.") do
      raise OptionParser::ParseError, "Only one of -T, -A and -Y are allowed." unless trim_method.nil?
      trim_method = [:trim_both]
    end
  end.parse!
  trim_method ||= [:trim_poly_a, :trim_poly_t] #by default trim both sides


  trimmed_count = [0,0,0]
  not_trimmed_count = [0,0,0]
  method_hash = {
    :trim_poly_a => 0,
    :trim_poly_t => 1,
    :trim_both=> 2,
  }

  Bio::FlatFile.open(ARGV[0]).each do |seq|
    trimmed = seq.seq
    trim_method.each do |method|
      last = trimmed
      trimmed = Bio::Trimpoly.send(method, last)

      # record how many have been trimmed
      index = method_hash[method]
      if last == trimmed
        not_trimmed_count[index] += 1
      else
        trimmed_count[index] += 1
      end
    end
    
    puts ">#{seq.definition}"
    puts trimmed
  end


  $stderr.puts "polyA: Trimmed #{trimmed_count[0]}, left #{not_trimmed_count[0]} untrimmed"
  $stderr.puts "polyT: Trimmed #{trimmed_count[1]}, left #{not_trimmed_count[1]} untrimmed"
  $stderr.puts "both: Trimmed #{trimmed_count[2]}, left #{not_trimmed_count[2]} untrimmed"
end