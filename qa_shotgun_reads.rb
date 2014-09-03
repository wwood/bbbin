#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'bio-commandeer'
require 'tempfile'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <arguments>

    preps illumina reads off the shotgun sequencer\n\n"

  opts.on("-1 FASTQ_GZ", "forward raw reads [required]") do |arg|
    options[:fwd] = arg
  end
  opts.on("-2 FASTQ_GZ", "forward raw reads [required]") do |arg|
    options[:rev] = arg
  end
  opts.on("-A ADAPTER", "forward adapter [required]") do |arg|
    options[:fwd_adapter] = arg
  end
  opts.on("-B ADAPTER", "reverse adapter [required]") do |arg|
    options[:rev_adapter] = arg
  end

  opts.on("-o BASENAME", "output basename [required]") do |arg|
    options[:output_basename] = arg
  end


  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end; o.parse!
if ARGV.length != 0 or options.length != 7
  $stderr.puts options
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME); log.outputters[0].formatter = Log4r::PatternFormatter.new(:pattern => "%5l %c %d: %m", :date_pattern => '%d/%m %T')

Tempfile.open(['fwd_tmp','.fq.gz']) do |fwd_tmp|
  Tempfile.open(['rev_tmp','.fq.gz']) do |rev_tmp|
    Bio::Commandeer.run "SeqPrep -f #{options[:fwd]} -r #{options[:rev]} -1 #{fwd_tmp.path} -2 #{rev_tmp.path} -A #{options[:fwd_adapter]} -B #{options[:rev_adapter]} -s #{options[:output_basename]}.seqprep_merged.fastq.gz", :log => true
    Bio::Commandeer.run "nesoni clip --quality 20 --homopolymers yes --length 30 #{options[:output_basename]} pairs: #{fwd_tmp.path} #{rev_tmp.path}", :log => true
  end
end
