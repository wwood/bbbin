#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'csv'



SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = { #TODO set defaults here
  :logger => 'stderr',
}
o = OptionParser.new do |opts|
  #TODO Fill in usage, description and option parsing below
  opts.banner = "
    Usage: #{SCRIPT_NAME} <keyword>
    
    Output a list of NCBI sequence read archive (SRA) accession numbers where a certain word is found in the study's abtract\n"

  
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
if ARGV.length != 1 #TODO require a set number of arguments?
  $stderr.puts o
  exit 1
end
keyword = ARGV[0]

# Setup logging
Bio::Log::CLI.logger(options[:logger]) #bio-logger defaults to STDERR not STDOUT, I disagree
log = Bio::Log::LoggerPlus.new(LOG_NAME)
Bio::Log::CLI.configure(LOG_NAME)


Bio::SRA::Tables::SRA.where(
  'study_abstract like ? or study_abstract like ?','%16S%','%16s%'
  ).where(:platform => 'LS454').where(
  'study_abstract like ? or study_accession like ?',"%#{keyword}%","%{keyword.capitalise}"
  ).where(:library_strategy => 'AMPLICON').all.each do |r|
  
  puts r.run_accession
end 
