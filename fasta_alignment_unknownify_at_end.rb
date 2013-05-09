#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'bio'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <aligned_fasta_file>

    Take an input aligned FASTA format file, and replace any dashes (which represent gaps) at the start and ends of each sequence with dots (representing missing)\n\n"


  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end; o.parse!
if ![0,1].include?(ARGV.length)
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)

Bio::FlatFile.foreach(Bio::FastaFormat, ARGF) do |seq|
  my_seq = seq.seq

  starting_dashes = my_seq.match(/^(\-+)(.+)/)
  unless starting_dashes.nil?
    log.debug "found #{starting_dashes[1]}" if log.debug?
    my_seq = '.'*starting_dashes[1].length + starting_dashes[2]
  end

  ending_dashes = my_seq.match(/(.+?)(-+)$/)
  unless ending_dashes.nil?
    log.debug "found #{ending_dashes}" if log.debug?
    my_seq = ending_dashes[1] + '.'*ending_dashes[2].length
  end

  puts '>'+seq.definition
  puts my_seq
end
