#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
$LOAD_PATH.unshift(File.join(ENV['HOME'],'git','bioruby-sra','lib'))
require 'bio-sra'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
}
o = OptionParser.new do |opts|
  #TODO Fill in usage, description and option parsing below
  opts.banner = "
    Usage: #{SCRIPT_NAME} <accession1> [<accession2> ..]
    
    output a list of accession numbers associated with one or more input accession numbers, so you can use google scholar etc.\n"
  # Example option
  # opts.on("-d", "--db", "description [default: #{options[:eg]}]") do |f|
    # options[:operation] = OVERALL
  # end
  
  # logger options
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") do |q|
    Bio::Log::CLI.trace('error')
  end
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") do | name |
    options[:logger] = name
  end
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG") do | s |
    Bio::Log::CLI.trace(s)
  end
end
o.parse!
if ARGV.length == 0 #TODO require a set number of arguments?
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]) #bio-logger defaults to STDERR not STDOUT, I disagree
log = Bio::Log::LoggerPlus.new(LOG_NAME)
Bio::Log::CLI.configure(LOG_NAME)


# Connect to the database
Bio::SRA.connect

# get a list of study, sample, etc ids
# uniq them
# print them all out joined by ' OR ' 
runs = ARGV.collect do |acc|
  Bio::SRA::Tables::SRA.accession(acc).all.collect do |sra|
    [sra.submission_accession, sra.study_accession, sra.sample_accession, sra.experiment_accession, sra.run_accession]
  end
end
puts runs.flatten.uniq.join(' OR ')
