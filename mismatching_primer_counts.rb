#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'bio-ipcress'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
  :forward => 'AAACTYAAAKGAATTGRCGG',
  :reverse => 'ACGGGCGGTGWGTRC',
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <arguments>

    Description of what this program does...\n\n"

  opts.separator "\nRequired:\n\n"
  opts.on("-i", "--input FASTA_FILE", "fasta file to search [required]") do |arg|
    options[:input] = arg
  end

  opts.separator "\nOptional:\n\n"
  opts.on("-f", "--forward PRIMER", "forward primer [default: #{options[:forward]}]") do |arg|
    options[:forward] = arg
  end
  opts.on("-r", "--input PRIMER", "forward primer [default: #{options[:reverse]}]") do |arg|
    options[:reverse] = arg
  end

  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end; o.parse!
if ARGV.length != 0
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)


primer_set = Bio::Ipcress::PrimerSet.new(
  options[:forward], options[:reverse]
  )

results = Bio::Ipcress.run(
  primer_set,
  options[:input], :mismatches => 1)

results.select! do |res|
  res.recalculate_mismatches_from_alignments == [0,0]
end

log.info "Found #{results.length} hits total"




