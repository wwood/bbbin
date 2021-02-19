#!/usr/bin/env ruby

require 'optparse'
# require 'bio-logger'
require 'bio'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :whole_fasta => false,
  :logger => 'stderr',
  :log_level => 'info',
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <fasta_file>

Calculate the GC content of each contig in the fasta file\n\n"


  opts.on("-w", "--whole", "Return the GC of the entire input file [default: report per-sequence]") {
    options[:whole_fasta] = true
  }
  opts.on('--interval NUM', Integer, "Return GC of each interval in a sequence e.g. every 1000 bp") do |d|
    options[:interval] = d
  end

  # logger options
  # opts.separator "\nVerbosity:\n\n"
  # opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  # opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  # opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end; o.parse!
if ARGV.length > 1
  $stderr.puts o
  exit 1
end
# Setup logging
# Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME); log.outputters[0].formatter = Log4r::PatternFormatter.new(:pattern => "%5l %c %d: %m", :date_pattern => '%d/%m %T')

raise "cant be whoel seq and interval" if options[:whole_fasta] and options[:interval]

total_seq = ''
Bio::FlatFile.foreach(ARGF) do |seq|
  if options[:whole_fasta]
    total_seq += seq.naseq.to_s
  elsif options[:interval]
    current = 0
    interval = options[:interval]
    while current + interval <= seq.naseq.length
      puts [
        seq.definition.split(/\s/)[0],
        current+1,
        current+interval,
        seq.naseq.subseq(current+1, current+interval).gc_content.to_f
      ].join("\t")
      current += interval
    end
  else
    puts [
      seq.definition.split(/\s/)[0],
      seq.naseq.gc_content.to_f,
    ].join("\t")
  end
end

if options[:whole_fasta]
  puts Bio::Sequence::NA.new(total_seq).gc_content.to_f
end
