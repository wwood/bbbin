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

    Given a distance matrix output from ExpressBetaDiversity, plot
    a grid of their relationships\n\n"

  opts.on("-i", "--input PATH", "path to ebd distance matrix file [required]") do |arg|
    options[:input_path] = arg
  end
  opts.on("--order FILE", "File containinng newline separated sample names in the order of the output [default: order from distance matrix]") do |arg|
    options[:order_file] = arg
  end
#  opts.on("-o", "--output-path PATH", "path to output file [required]") do |arg|
#    options[:output_path] = arg
#  end

  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end; o.parse!
if ARGV.length != 0 or options[:input_path].nil?
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)

# Read in the ebd matrix
ebd_matrix = Bio::EBD::DistanceMatrix.parse_from_file(options[:input_path])
log.info "Read in distance matrix containing #{ebd_matrix.number_of_samples} samples"

order = ebd_matrix.sample_names
if options[:order_file]
  order = File.open(options[:order_file]).read.split("\n").collect{|r| r.strip}
  log.info "Read in #{order.length} samples for ordering / subsetting"
end


puts '['
order.each_with_index do |name1,i|
  print '['
  order.each_with_index do |name2,j|
    if i==j
      print '0'
    else
      print 1-ebd_matrix.distance(name1,name2)
    end
    print ','
  end
  puts '],'
end
puts ']'
