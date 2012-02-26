#!/usr/bin/env ruby

# Apply a sliding window over a fasta file, outputing the number of glycines (and some other characteristics) found maximally in those proteins.

require 'rubygems'
require 'bio'

module Bio
  class SequenceWindowDescriptor
    attr_reader :maximum_counts, :maximum_sequences
    def calculate(sequence_object, window_size)
      # initialise maximums
      @maximum_counts = {}
      @maximum_sequences = {}
      [:gly].each do |sym|
        @maximum_counts[sym] = 0
        @maximum_sequences[sym] = ''
      end

      sequence_object.window_search(window_size) do |str|
        num = 0
        str.scan(/g/i) {num += 1}
        if num >= @maximum_counts[:gly]
          @maximum_counts[:gly] = num
          @maximum_sequences[:gly] = str
        end
      end
    end
  end
end

if __FILE__ == $0
  require 'optparse'

  # Parse cmd line options
  USAGE = "Usage: gly_sliding_window.rb <fasta_filename>"
  options = {
    :window_size => 25,
  }
  o = OptionParser.new do |opts|
    opts.banner = USAGE

    opts.on("-w", "--window-size SIZE", "Length of the window to be used") do |v|
      window = v.to_i
      unless window > 0
        raise Exception, "Unexpected window size specified: #{v} - it must be greater than 0 residues long!"
      end
      options[:window_size] = window
    end
  end
  o.parse!

  puts %w(name window_size max_number max_sequence).join("\t")
  Bio::FlatFile.foreach(ARGF) do |seq|
    desc = Bio::SequenceWindowDescriptor.new
    desc.calculate(seq.seq, options[:window_size])
    puts [seq.definition, options[:window_size], desc.maximum_counts[:gly], desc.maximum_sequences[:gly]].join("\t")
  end
end