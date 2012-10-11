#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'progressbar'

if __FILE__ == $0 #needs to be removed if this script is distributed as part of a rubygem
  SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')
  
  # Parse command line options into the options hash
  options = {
    :logger => 'stderr',
  }
  o = OptionParser.new do |opts|
    opts.banner = "
      Usage: #{SCRIPT_NAME} -p pileup_file -f identifiers_file
      
      Takes a pileup file and only keeps a subset of reference sequences, to reduce the size.\n\n"
      
    opts.on("-p", "--pileup PILEUP_FILE", "Pileup file to be subsetted [required]") do |arg|
      options[:pileup_file] = arg
    end
    opts.on("-f", "--references-to-keep IDENTIFIERS_FILE", "File ofreference sequence IDs that are to be kept [required]") do |arg|
      options[:identifiers_file] = arg
    end
    
    # logger options
    opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {Bio::Log::CLI.trace('error')}
    opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
    opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| Bio::Log::CLI.trace(s)}
  end
  o.parse!
  if ARGV.length != 0 or options[:pileup_file].nil? or options[:identifiers_file].nil?
    $stderr.puts o
    exit 1
  end
  # Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
  Bio::Log::CLI.logger(options[:logger]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)
  
  
  # read in the list of identifiers
  hash = {}
  ids = File.open(options[:identifiers_file]).each_line do |line|
    r = line.strip
    unless r.length == 0
      log.warn "reference '#{r}' specified multiple times in the identifiers file, I'll only print it once" if hash.key?(r)
      hash[r] = true
    end
  end
  log.info "Cached #{hash.length} reference identifiers from #{options[:identifiers_file]}"
  
  # Print out the subset
  total = `wc -l '#{options[:pileup_file]}'`.to_i
  progress = ProgressBar.new('pileup_subset', total)
  File.open(options[:pileup_file]).each_line do |line|
    row = line.split("\t")
    
    if row.length > 0 and hash.key?(row[0])
      puts line
    end
    progress.inc
  end
  progress.finish
  
end #end if running as a script