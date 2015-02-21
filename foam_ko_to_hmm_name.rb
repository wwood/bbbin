#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <arguments>

Takes a multi-hmm FOAM file, and outputs a list of KEGG KO (e.g. KO:K10759) to HMM name (e.g. HMMsoil99735) on stdout.\n\n"

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


# HMMER3/b [3.0 | March 2010]
# NAME  KO:K10759 3.1.1.20
# TC    557.6 556.8;
# ACC   HMMsoil99735
# DESC  KO:K10759 3.1.1.20
# LENG  453
# ALPH  amino

log.info "Parsing.."


last_name = nil
num_pairs = 0
ARGF.each_line do |line|
  if line[0..2] == 'ACC'
    last_name = line.split(/ +/)[1]
  elsif line[0..3] == 'DESC'
    line.strip.split(/ +/).each_with_index do |e, i|
      next if i==0
      if e.match(/^KO:/)
        puts [
          last_name,
          e
          ].join("\t")
        num_pairs += 1
      end
    end
  end
end
log.info "Printed #{num_pairs} pairs"
