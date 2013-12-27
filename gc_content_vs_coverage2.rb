#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
#require 'bio-samtools'
require 'bio'
require 'progressbar'
require 'csv'
require 'bio-commandeer'
require 'hopcsv'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :window_length => 100,
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <arguments>

    Take an (indexed) bam file and and plot the coverages of windows from the base sequence, recording the GC content of the window\n\n"

  opts.separator "\nRequired:\n\n"
  opts.on('-f',"--fasta FILENAME", "Fasta file of the reference. Must be indexed [required]") do |arg|
    options[:fasta_file] = arg
  end
  opts.on('-b',"--bam FILENAME", "BAM file to get input from. Must be indexed [required]") do |arg|
    options[:bam_file] = arg
  end
  opts.on("--depths FILENAME", "Output from samtools depth [required]") do |arg|
    options[:depth_file] = arg
  end

  opts.separator "\nOptional:\n\n"
  opts.on("--window-size LENGTH", "Length of window to cover [default: #{options[:window_length]}]") do |arg|
    options[:window_length] = arg.to_i
    raise unless options[:window_length]>0
  end

  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {Bio::Log::CLI.trace('error')}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| Bio::Log::CLI.trace(s)}
end; o.parse!
if ARGV.length != 0 or options[:fasta_file].nil or (optoins[:bam_file].nil? and options[:depth_file].nil?)
  $stderr.puts o
  exit 1
end
# Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
Bio::Log::CLI.logger(options[:logger]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)


# For each window of the reference sequence
# sam = Bio::DB::Sam.new(
  # :bam => options[:bam_file],
  # :fasta => options[:fasta_file],
# )
# sam.load_index
# sam.load_reference
# sam.open
# log.debug "Loaded SAM object: #{sam.to_s}"

puts %w(
  sequence
  start
  stop
  gc
  coverage
  positions_covered
  ).join "\t"
total_zero_coverage_points = 0

# Cache sequences
sequences = {}
Bio::FlatFile.foreach(options[:fasta_file]) do |s|
  identifier = s.definition.split(/\s/)[0]
  sequences[identifier] = s.seq.seq
end

prev_coverages = []
first_seq_name = nil
i = 1
next_ceiling = options[:window_length]
floor_position = i

progress = ProgressBar.new('1',sequences.values[0].length/options[:window_length])

depth_file = options[:depth_file]
tempfile = nil
if depth_file.nil?
  tempfile = Tempfile.new('depth')
  cmd = "samtools depth #{options[:bam_file].inspect} >#{tempfile.path}"
  Bio::Commandeer.run cmd, :log => true
  depth_file = tempfile.path
end

HopCSV.foreach(depth_file, "\t") do |row|
  seq_name = row[0]
  if first_seq_name.nil?
    first_seq_name = row[0]
  elsif first_seq_name != seq_name
    raise "Cannot handle multiple reference sequences in the bam file at the moment"
  end

  pos = row[1].to_i
  cov = row[2].to_i

  while pos > next_ceiling
    progress.inc

    # Fill in the zeroes if they occur
    num_filled = 0
    while prev_coverages.length < options[:window_length]
      prev_coverages.push 0
      num_filled += 1
    end
    log.debug "Filled in #{num_filled} 0 coverage points in this window" if log.debug? and num_filled > 0
    total_zero_coverage_points += num_filled

    puts [
      seq_name,
      floor_position,
      next_ceiling,
      prev_coverages.reduce(:+).to_f/prev_coverages.length,
      Bio::Sequence::NA.new(sequences[seq_name][floor_position-1,options[:window_length]]).gc_percent,
      prev_coverages.length,
    ].join("\t")

    log.error "PROGRAMMING ERROR Coverage length between #{floor_position} and #{next_ceiling} was #{prev_coverages.length}" unless prev_coverages.length == options[:window_length]

    floor_position = next_ceiling+1
    prev_coverages = []
    next_ceiling = floor_position+options[:window_length]-1

  end

  i+=1
  prev_coverages.push cov
end
progress.finish
$stderr.puts "Found #{total_zero_coverage_points} places in the genome that had zero coverage"




