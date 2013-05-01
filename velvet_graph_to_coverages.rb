#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'bio-velvet'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <velvet_graph_file>

    Takes a velvet graph file and outputs a tab separated file of each node's coverage in that graph\n\n"

  #opts.on("-e", "--eg ARG", "description [default: #{options[:eg]}]") do |arg|
  #  options[:example] = arg
  #end

  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end; o.parse!
if ARGV.length != 1
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)


#TODO what are you looking at me for? This is your script. Do something.
graph_file = ARGV[0]
graph = Bio::Velvet::Graph.parse_from_file graph_file

graph.nodes.each do |node|
  puts [
    node.node_id,
    node.coverages
  ].flatten.join "\t"
end
