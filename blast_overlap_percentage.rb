#!/usr/bin/env ruby
require 'csv'
require 'optparse'
require 'bio-logger'
require 'bio-genomic-interval'

if __FILE__ == $0 #needs to be removed if this script is distributed as part of a rubygem
  SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')
  
  # Parse command line options into the options hash
  options = {
    :logger => 'stderr',
  }
  o = OptionParser.new do |opts|
    opts.banner = "
      Usage: #{SCRIPT_NAME} -b blast_csv_file
      
      Takes a blast tab-separated values file, and determines the length of query sequence that is overlapped by at least one hit sequence\n\n"
      
    opts.on("-b", "--blast <blast_file>", "blast file to work off [required]") do |arg|
      options[:blast_file] = arg
    end
    opts.on("-s", "--subject", "Compute for the subject, not the query sequence [default: false (use query, not subject)]") do
      options[:use_subject] = true
    end

    # logger options
    opts.separator "\nVerbosity:\n\n"
    opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {Bio::Log::CLI.trace('error')}
    opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
    opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| Bio::Log::CLI.trace(s)}
  end; o.parse!
  if ARGV.length != 0 or options[:blast_file].nil?
    $stderr.puts o
    exit 1
  end
  # Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
  Bio::Log::CLI.logger(options[:logger]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)
  
  
  intervals = []
  CSV.foreach(options[:blast_file], :col_sep => "\t") do |row|
    if options[:use_subject]
      intervals.push Bio::GenomicInterval.parse("#{row[1]}:#{row[8]}-#{row[9]}")
    else
      intervals.push Bio::GenomicInterval.parse("#{row[0]}:#{row[6]}-#{row[7]}")
    end
  end
  
  # Merge overlaps
  overlap_just_removed = true
  while overlap_just_removed
    overlap_just_removed = false
    intervals.combination(2).each do |combo|
      interval1 = combo[0]
      interval2 = combo[1]
      
      if interval1.overlapped?(interval2)
        intervals = [
          intervals.reject{|interval| interval==interval1 or interval==interval2},
          interval1.expand(interval2)          
        ].flatten
        overlap_just_removed = true
        log.debug "#{interval1} overlaps with #{interval2}"
        break
      end
    end
  end
  
  # Print total length of overlap
  puts intervals.collect{|interval| interval.length}.reduce(:+)
end #end if running as a script