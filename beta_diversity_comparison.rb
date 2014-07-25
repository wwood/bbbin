#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'bio-express_beta_diversity'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <arguments>

    Output distances in an Express Beta Diversity .diss file between two sets of sample IDs\n\n"

  opts.on("--ebd DISS_FILE", "Express Beta Diversity .diss file [required]") do |arg|
    options[:diss] = arg
  end
  opts.on("--sample-ids1 IDs", Array, "First comma-separated list of sample IDs for comparison [required]") do |arg|
    options[:sample_ids1] = arg
  end
  opts.on("--sample-ids2 IDs", Array, "Second comma-separated list of sample IDs for comparison [required]") do |arg|
    options[:sample_ids2] = arg
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


log.info "Reading in EBD distance file `#{options[:diss]}' .."
ebd = Bio::EBD::DistanceMatrix.parse_from_file options[:diss]
log.info "Read in distances between #{ebd.sample_names.length} samples"

options[:sample_ids1].each do |s1|
  raise "Sample #{s1} not found in EBD file" unless ebd.sample_names.include?(s1)
  options[:sample_ids2].each do |s2|
    next if s1==s2
    raise "Sample #{s2} not found in EBD file" unless ebd.sample_names.include?(s2)
    puts [
      s1,
      s2,
      ebd.distance(s1, s2)
    ].join("\t")
  end
end
