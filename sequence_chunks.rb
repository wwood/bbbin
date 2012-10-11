#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'bio'

if __FILE__ == $0 #needs to be removed if this script is distributed as part of a rubygem
  SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')
  
  # Parse command line options into the options hash
  options = {
    :logger => 'stderr',
    :chunk_size => 500,
    :step_size => 250,
  }
  o = OptionParser.new do |opts|
    opts.banner = "
      Usage: #{SCRIPT_NAME} <fasta_file>
      
      Takes an input fasta file, and separates each sequence it into overlapping pieces. Prints the chunks on stdout.\n\n"
      
    opts.on("-c", "--chunk-size ARG", "Size of the chunks to output [default: #{options[:chunk_size]}]") do |arg|
      options[:chunk_size] = arg.to_i
      raise unless options[:chunk_size] > 0
    end
    opts.on("-s", "--step-size ARG", "Size of the chunks to output [default: #{options[:step_size]}]") do |arg|
      options[:step_size] = arg.to_i
      raise unless options[:step_size] > 0
    end
    
    # logger options
    opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {Bio::Log::CLI.trace('error')}
    opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
    opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| Bio::Log::CLI.trace(s)}
  end
  o.parse!
  if ![0,1].include?(ARGV.length)
    $stderr.puts o
    exit 1
  end
  # Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
  Bio::Log::CLI.logger(options[:logger]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)
  
  window_size = options[:chunk_size]
  step_size = options[:step_size]
    
  print_seq = lambda do |fasta, subseq, start, number_in_sequence, total_fragment_count|
    first_name = "#{fasta.definition.split(/\s/)[0]}_#{number_in_sequence}of#{total_fragment_count}_#{start+1}_#{start+window_size}"
    print '>'
    puts first_name
    puts subseq
  end
  
  Bio::FlatFile.foreach(ARGF) do |seq|
    fragment_number = 1
    total_fragment_count = seq.length.to_f/step_size

    is_exact_multiple_of_length = ((seq.length.to_f/step_size).round == seq.length/step_size)
    total_fragment_count += 1 unless is_exact_multiple_of_length #add one for leftovers unless it is exactly on the money 

    total_fragment_count = total_fragment_count.to_i
    
    last_step = 0
    0.step(seq.length - window_size, step_size) do |i| 
      subseq = seq.seq[i, window_size]
      print_seq.call(seq, subseq, i, fragment_number, total_fragment_count)
      last_step = i
      fragment_number += 1
    end
    start = seq.length-window_size
    unless is_exact_multiple_of_length
      leftover = seq.seq[start .. -1]
      print_seq.call(seq, leftover, start, fragment_number, total_fragment_count)
    end
  end
end #end if running as a script