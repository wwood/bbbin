#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'bio-commandeer'
require 'bio'
require 'pry'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <arguments>

    Description of what this program does...\n\n"

  opts.on("--protein-file ARG", "description [required]") do |arg|
    options[:protein_file] = arg
  end
  opts.on("--nucleotides-file ARG", "description [required]") do |arg|
    options[:nucleotide_file] = arg
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
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME); log.outputters[0].formatter = Log4r::PatternFormatter.new(:pattern => "%5l %c %d: %m", :date_pattern => '%d/%m %T')




cmd = "tblastn -query #{options[:protein_file]} -subject #{options[:nucleotide_file]} -outfmt 6"

out = Bio::Commandeer.run(cmd)
if out.split("\n").length < 1
  raise
else
  splits = out.split("\n")[0].split "\t"
  #Pfer_Contig00027_12027_6_213	Pfer_Contig00027	100.00	241	0	0	1	241	12749	12027	2e-156	  480
  start = splits[8].to_i
  stop = splits[9].to_i
  if start > stop
    s = stop
    stop = start
    start = s
  end
  name = splits[1]

  names_to_sequences = {}
  Bio::FlatFile.open(File.open(options[:nucleotide_file])).entries.each do |s|
    names_to_sequences[s.definition.split(/\s/)[0]] = s.seq.seq
  end

  start = start - 500
  start = 1 if start < 1
  stop = stop + 500
  seq = names_to_sequences[name]
  raise if seq.nil?
  stop = seq.length if stop > seq.length

  puts ">#{name}_#{start}:#{stop}"
  puts seq[(start-1)..stop]
end
