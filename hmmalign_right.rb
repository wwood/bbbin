#!/usr/bin/env ruby
require 'pp'
require 'optparse'
require 'bio-logger'
require 'bio'
require 'tempfile'
require 'bio-stockholm'
require 'pry'
SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <FASTA_FILE>

    Uses hmmalign, but removes unaligned columns, which often aren't helpful but instead just get in the way\n\n"

  opts.on("--hmm HMMFILE", "HMM file to align to [required]") do |arg|
    options[:hmm_file] = arg
  end
  opts.on("--fasta FASTA_FILE", "Fasta file of sequences to align [required]") do |arg|
    options[:fasta] = arg
  end

  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end; o.parse!
if ARGV.length > 0 or options[:hmm_file].nil? or options[:fasta].nil?
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)



Tempfile.open(['hmmalign_right_out','.sto']) do |out_stockholm|
  cmd = "hmmalign --trim #{options[:hmm_file]} #{options[:fasta]} >#{out_stockholm.path}"
  log.debug "Running: #{cmd}" if log.debug?
  print `#{cmd}`

  stocks = Bio::Stockholm::Reader.parse_from_file out_stockholm.path
  stocks[0].records.each do |identifier, record|
    puts ">#{identifier}"
    puts record.sequence.gsub(/[\.a-z]/,'')
  end
end
