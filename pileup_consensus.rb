#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'bio-pileup_iterator'


if __FILE__ == $0 #needs to be removed if this script is distributed as part of a rubygem
  SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')
  
  # Parse command line options into the options hash
  options = {
    :logger => 'stderr',
  }
  o = OptionParser.new do |opts|
    opts.banner = "
      Usage: #{SCRIPT_NAME} <pileup_file>
      
      From a pileup file, print the most prevalent base at each position.\n\n"

    # logger options
    opts.separator "\n\tVerbosity:\n\n"
    opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {Bio::Log::CLI.trace('error')}
    opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
    opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| Bio::Log::CLI.trace(s)}
  end; o.parse!
  if ARGV.length != 0
    $stderr.puts o
    exit 1
  end
  # Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
  Bio::Log::CLI.logger(options[:logger]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)
  
  
  Bio::DB::PileupIterator.new(ARGF).each do |pileup|
    # Pileup#consensus returns the consensus (most frequent) base from the pileup, 
    # if there are equally represented bases returns a string containing all equally 
    # represented bases in alphabetical order.
    #
    # Here I'll just pick a random one.
    consensus = pileup.consensus
    sample_index = (0...consensus.length).to_a.sample
    print consensus[sample_index]
  end
  puts
end #end if running as a script