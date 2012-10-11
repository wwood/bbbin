#!/usr/bin/env ruby

# Apply a sliding window over a fasta file, outputing the number of glycines (and some other characteristics) found maximally in those proteins.

require 'rubygems'
require 'bio'
require 'bio-signalp'
require 'csv'

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
    :signalp_file => nil,
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
    
    opts.on("-s", "--signalp_overview SIGNALP_SUMMRY_FILE", "Path to a bio-signalp derived summary file. Implies that only signalp-positive proteins are considered, and the gly sliding window cannot overlap with the signal peptide") do |v|
      options[:signalp_file] = v
    end
  end
  o.parse!
  
  signalp_cleavages = {}
  if options[:signalp_file]
    headering = true
    File.open(options[:signalp_file]).each_line do |line|
      if headering
        headering = false
        next
      end
      next if line.strip.length == 0
      # Name NN Prediction HMM Prediction  Predicted?  Cleavege site (if predicted)
      row = line.strip.split("\t")
      name = row[0]
      cleavage = row[4].to_i

      # Check for duplicates
      if signalp_cleavages[name]
        raise Exception, "Unexpectedly found sequence with name #{name.inspect} (at least) twice in the signalp_overview file, giving up"
      end
      signalp_cleavages[name] = cleavage
    end
    $stderr.puts "Recorded signalp predictions for #{signalp_cleavages.length} sequences, e.g. #{signalp_cleavages.to_a[0].inspect}"
  end

  headers = %w(name window_size max_number max_sequence)
  headers.push 'signalp?' if options[:signalp_file]
  puts headers.join("\t")
  
  Bio::FlatFile.foreach(ARGF) do |seq|
    desc = Bio::SequenceWindowDescriptor.new
    
    # By default use the whole sequence
    sequence = seq.seq
    name = seq.entry_id
    
    # If a signal peptide was predicted, cleave
    result = Bio::SignalP::Version3::Result.new
    if options[:signalp_file]
      # bit of a hack - artificially set the signalp value to true
      signalp_result = signalp_cleavages[name]
      if signalp_result.nil?
        $stderr.puts "Didn't find sequence #{name} in the signalp overview file, is that right?"
        signalp_result = 0
      end
      result.nn_D_prediction = (signalp_result != 0)
      result.nn_Ymax_position = signalp_result
      sequence = result.cleave(seq.seq)
    end
    
    desc.calculate(sequence, options[:window_size])
    print [seq.definition, options[:window_size], desc.maximum_counts[:gly], desc.maximum_sequences[:gly]].join("\t")
    if options[:signalp_file]
      print "\t"
      print signalp_cleavages[name] != 0
    end
    puts
  end
end
