#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'bio'
require 'set'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
  :min_piece_length => 1000,
  :chunk_removal_min => 0,
  :chunk_removal_max => 0,
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} [arguments] <fasta_file>

    Take a fasta file, and randomly split some number of them into two parts\n\n"

  opts.on("-n", "--number-to-split NUM", Integer, "split these number of sequences up [required]") do |arg|
    options[:number_to_split] = arg.to_i
  end
  opts.on("--min-length LENGTH", Integer, "don't split contigs up into pieces less than this long [default #{options[:min_piece_length]}]") do |arg|
    options[:min_piece_length] = arg.to_i
  end
  opts.on("--remove-chunks LENGTH1,LENGTH2", Array, "As well as breaking the sequence, remove some of it too. Choose a random amount of sequence between the two lengths given to remove. [default: don't remove]") do |args|
    if args.length != 2
      raise "Bad argument for --remove-chunks: #{args.inspect}"
    end
    options[:chunk_removal_min] = args[0].to_i
    options[:chunk_removal_max] = args[1].to_i
    unless options[:chunk_removal_min] <= options[:chunk_removal_max]
      raise "max chunk removal size smaller than min, I'm confused!"
    end
  end

  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end; o.parse!
if ARGV.length != 1
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)

fasta_file = ARGV[0]

# Count the number of sequences in the fasta file that are greater than twice the minimum length
sufficient_length = 2*options[:min_piece_length]+1+options[:chunk_removal_max]
num_seqs_with_sufficient_length = 0
log.debug "Reading from fasta file #{fasta_file}"
Bio::FlatFile.foreach(fasta_file) do |seq|
  log.debug "Considering a sequence #{seq.definition}, length #{seq.seq.length}" if log.debug?
  if seq.seq.length >= sufficient_length
    num_seqs_with_sufficient_length += 1
    log.debug "Sequence accepted" if log.debug?
  end
end
log.info "Found #{num_seqs_with_sufficient_length} sequences of sufficient length"


# Ensure that the number to split is less than the number of entries
if num_seqs_with_sufficient_length < options[:number_to_split]
  raise "Insufficient number of sequences #{num_seqs_with_sufficient_length} to take #{options[:number_to_split]}"
end

# Choose the fasta file entries that are to be split
chosen_indices = Set.new (0...num_seqs_with_sufficient_length).to_a.shuffle[0...options[:number_to_split]]

# Iterate through the fasta file again, splitting reads as necessary
random_gen = Random.new
sufficient_index = 0;
Bio::FlatFile.foreach(fasta_file) do |seq|
  # If this read is greater than the min length
  if seq.seq.length >= sufficient_length
    num_seqs_with_sufficient_length =+ 1

    # If this is chosen to be split, split it
    if chosen_indices.include?(sufficient_index)
      chunk_removal_size = random_gen.rand(options[:chunk_removal_max]-options[:chunk_removal_min]) + options[:chunk_removal_min]
      split_point = random_gen.rand(seq.seq.length-2*options[:min_piece_length]-chunk_removal_size) + options[:min_piece_length]
      seq1 = seq.seq[0...split_point]
      seq2 = seq.seq[(split_point+chunk_removal_size)...seq.seq.length]
      puts ">#{seq.definition}_break1_#{split_point}_#{chunk_removal_size}"
      puts seq1
      puts ">#{seq.definition}_break2_#{split_point}_#{chunk_removal_size}"
      puts seq2
    else
      # Else just print it
      puts seq.to_s
    end

    # Increment the counter of sufficient length contigs
    sufficient_index += 1
  else
    puts seq.to_s
  end

end
