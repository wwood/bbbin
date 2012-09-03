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
      Usage: #{SCRIPT_NAME} <genes1.nucleotide.fasta> <genes2.nucleotide.fasta>
      
      Takes two sets of gene predictions, and outputs a table that gives the correspondences between them. Matching is based purely on sequence.\n\n"

    # logger options
    opts.separator "\nVerbosity:\n\n"
    opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {Bio::Log::CLI.trace('error')}
    opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
    opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| Bio::Log::CLI.trace(s)}
  end; o.parse!
  if ARGV.length != 2
    $stderr.puts o
    exit 1
  end
  # Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
  Bio::Log::CLI.logger(options[:logger]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)
  
  
  genes1_sequence_to_id_hash = {}
  Bio::FlatFile.foreach(ARGV[0]) do |s|
    ident = s.definition
    seq = s.seq
    
    if genes1_sequence_to_id_hash[seq]
      log.warn "Sequence clash found in the first set, between '#{genes1_sequence_to_id_hash[seq]}' and '#{ident}', ignoring the second one"
    else
      genes1_sequence_to_id_hash[seq] = ident
    end
  end
  
  genes2_sequence_to_id_hash = {}
  genes1_matched = []
  Bio::FlatFile.foreach(ARGV[1]) do |s|
    ident = s.definition
    seq = s.seq
    
    if genes2_sequence_to_id_hash[seq]
      log.warn "Sequence clash found in the second set, between '#{genes2_sequence_to_id_hash[seq]}' and '#{ident}', ignoring the second one"
    else
      if genes1_sequence_to_id_hash[seq]
        puts [
          genes1_sequence_to_id_hash[seq],
          ident,
        ].join("\t")
        genes2_sequence_to_id_hash[seq] = ident
        genes1_matched.push genes1_sequence_to_id_hash[seq]
      else
        log.warn "Unable to match from 2nd fasta #{ident}"
      end
    end
  end
  
  # Output all those from the first file unmatched
  genes1_sequence_to_id_hash.values.each do |ident|
    unless genes1_matched.include?(ident)
      log.warn "Unable to match from 1st fasta #{ident}"
    end
  end
end #end if running as a script