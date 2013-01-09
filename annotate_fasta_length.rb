#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'bio'

if __FILE__ == $0 #needs to be removed if this script is distributed as part of a rubygem
  SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')
  
  # Parse command line options into the options hash
  options = {
    :logger => 'stderr',
  }
  o = OptionParser.new do |opts|
    opts.banner = "
      Usage: #{SCRIPT_NAME} <fasta_file>
      
      Takes a fasta file, and then adds the length of sequence to the identifier e.g. '>myseq' => '>myseq_len200346'.
      
      Outputs the modified fasta file to STDOUT\n\n"

    # logger options
    opts.separator "\nVerbosity:\n\n"
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
  
  
  Bio::FlatFile.foreach(ARGF) do |e|
    splits = e.definition.split(/\s/)
    extras = ''
    if splits.length > 1
      extras = splits[1...splits.length].join(' ')
    end
    length = e.seq.length
    first_name = splits[0]+'_len'+length.to_s
    if splits.length > 1
      e.definition = [
        first_name,
        extras
      ].join(' ')
    else
      e.definition = first_name
    end
    puts e
  end
  
end #end if running as a script