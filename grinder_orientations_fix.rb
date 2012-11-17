#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'bio-faster'
require 'bio'

if __FILE__ == $0 #needs to be removed if this script is distributed as part of a rubygem
  SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')
  
  # Parse command line options into the options hash
  options = {
    :logger => 'stderr',
  }
  o = OptionParser.new do |opts|
    opts.banner = "
      Usage: #{SCRIPT_NAME}
      
      Grinder, a software for simulating next-generation reads, has somewhat of a bug in it when it is generating forward pairs compared to reverse pairs.
      
      The forward strand pairs are inward facing, the reverse strand pairs are outward facing. This script changes the reverse pairs so that they are inward facing.
      
      Give a fastq file on stdin, and a fixed fastq will put printed to stdout.
      
      NOTE: Only the sequence and the quality is fixed, not the headers. The headers will continue to be wrong.
      \n\n"

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
  
  Bio::Faster.new(:stdin).each_record(:quality => :raw) do |sequence_header, sequence, quality|
    # If this is part of a pair simulated from the forward strand, just spit it straight back out.
    if (!sequence_header.match(/position=complement/) and sequence_header.match(/^\S+\/1 /)) or
      (sequence_header.match(/position=complement/) and sequence_header.match(/^\S+\/2 /))
      
      puts '@'+sequence_header
      puts sequence
      puts '+'
      puts quality

    else
      # This is simulated from the reverse strand. Reverse the reads and their orientations
      
      # Check for sanity
      if (!sequence_header.match(/position=complement/) and sequence_header.match(/^\S+\/2 /)) or
        (sequence_header.match(/position=complement/) and sequence_header.match(/^\S+\/1 /))
        
        puts '@'+sequence_header
        puts Bio::Sequence::NA.new(sequence).reverse_complement.to_s.upcase
        puts '+'
        puts quality.reverse
        
      else
        raise "Unexpected sequence header found: #{sequence_header}"
      end
    end
  end
end #end if running as a script