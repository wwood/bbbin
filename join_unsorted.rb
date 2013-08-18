#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
  :expect_duplicate_keys => false,
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} file1 file2

    GNU join confuses me, so I'm going to write something that doesn't require sorting, which is the usual case

    Either file1 or file2 (but not both) can be '-' i.e. read from STDIN not a file

    \n\n"

  opts.on("--expect-duplicate-keys", "Don't complain & skip if duplicate keys are found in the 2nd file [default: #{options[:expect_duplicate_keys]}]") do |arg|
    options[:expect_duplicate_keys] = arg
  end

  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end; o.parse!
if ARGV.length != 2
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)

separator = "\t"
join_field = 0 #careful changing this because it is hard-coded in some places

file1_hash = {}

# Gather file1 up
first_io = (ARGV[0] == '-' ? $stdin : File.open(ARGV[0])) #read first file from stdin or from a file?
first_io.each_line do |line|
  splits = line.chomp.split(separator)
  key = splits[join_field].strip
  if file1_hash[key] #assume keys are unique, advise and ignore subsequents otherwise
    log.warn "Found duplicate key in file 1: `#{key}', ignoring after the first time"
    next
  else
    file1_hash[key] = splits
  end
end

# Go through file 2
file1_foundeds = {}
num_splits_in_file1 = nil
second_io = (ARGV[1] == '-' ? $stdin : File.open(ARGV[1])) #read second file from stdin or from a file?
second_io.each_line do |line|
  splits = line.chomp.split(separator)

  if num_splits_in_file1.nil?
    num_splits_in_file1 = splits.length #assume all lines have the same number of splits, fail otherwise
  elsif num_splits_in_file1 != splits.length
    raise Exception, "Unexpected number of splits over '#{separator}' found in the first file. I give up so I don't cause you pain later"
  end

  key = splits[join_field].strip
  file2_data = splits[1..splits.length-1]
  if file1_hash[key]
    if file1_foundeds[key].nil? or options[:expect_duplicate_keys]
      puts [file1_hash[key],file2_data].flatten.join(separator)
      file1_foundeds[key] ||= 0
      file1_foundeds[key] += 1
    else
      log.warn "Found key `#{key}' more than once in the second file. Fail."
    end
  else
  # Do nothing -they aren't joined
#    # Print empty cells for the first file
#    (num_splits_in_file1).times {print separator}
#    # Print the data for the second file
#    puts file2_data.join(separator)
  end
end

# Print out all the remaining ones that aren't in file 2
#file1_hash.each do |key, value|
#  puts value.join(separator) unless file1_foundeds[key]
#end
