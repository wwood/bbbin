#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'bio-commandeer'
require 'set'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
  :continue => false,
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <arguments>

    Assemble read sets combinatorially e.g all individual, then all pairs, then all groups of 3 etc.\n\n"

  opts.on("-1 FILE", Array, "fq.gz read 1 files [required]") do |arg|
    options[:firsts] = arg
  end
  opts.on("-2 FILE", Array, "fq.gz read 2 files [required]") do |arg|
    options[:seconds] = arg
  end
  opts.on("-t INT", Integer, "number of threads") do |arg|
    options[:threads] = arg
  end
  opts.on("-o PATH", '--output', "Base output folder") do |arg|
    options[:output] = File.absolute_path arg
  end
  opts.on("-c COMBINATION_COUNTS", '--combination-counts', Array, "Only do combinations of this many read sets") do |arg|
    options[:combination_counts] = arg.collect{ |i| i.to_i }
  end
  opts.on("--continue", "Continue with more combinations, skip output directory creation") do
    options[:continue] = true
  end

  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end; o.parse!
if ARGV.length != 0 or not options[:firsts] or not options[:seconds] or not options[:output] or not options[:threads]
  require 'pry'; binding.pry
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)

raise if options[:firsts].length != options[:seconds].length

read_sets = []
names = Set.new

if not options[:continue]
  if File.exists? options[:output]
    raise Exception, "Output directory already exists!"
  else
    Dir.mkdir options[:output]
  end
end

get_name = lambda do |read_path|
  File.basename(read_path)
end

options[:firsts].each_with_index do |r1, i|
  name = get_name.call r1
  raise "duplicate read file basename" if names.include? name
  names.add name
  r2 = options[:seconds][i]
  read_sets << [File.absolute_path(r1), File.absolute_path(r2)]
end

run_assembly = lambda do |read_set_combination, output_folder|
  read1 = read_set_combination.collect{ |r| "'#{r[0]}'" }.join(',')
  read2 = read_set_combination.collect{ |r| "'#{r[1]}'" }.join(',')

  Bio::Commandeer.run "megahit -t #{options[:threads].to_s} -1 #{read1} -2 #{read2} -o '#{output_folder}' --tmp /tmp", :log => log
end

combination_counts = options[:combination_counts]
combination_counts ||= (1..read_sets.length)
combination_counts.each do |num|
  log.info("Running combinations of #{num}")

  base_output = File.join(options[:output], "combination#{num}")
  Dir.mkdir base_output
  read_sets.combination(num).each do |read_set_combination|
    output = File.join(base_output, read_set_combination.collect{ |r| get_name.call(r[0]) }.join('_with_'))
    run_assembly.call read_set_combination, output
  end
end
log.info "Finished running assemblies"


