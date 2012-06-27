#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'bio'


module Bio
  module Alignment
    class OriginalAlignment
      def trim_uninformative_columns
        newseqs = Bio::Alignment.new
        keys.each do |key|
          newseqs.add_seq('', key)
        end
        
        each_site do |site|
          alphabet_counts = {}
          site.each do |s|
            alphabet_counts[s] ||= 0
            alphabet_counts[s] += 1
          end
          
          keys = alphabet_counts.keys
          if alphabet_counts.length > 2 or (alphabet_counts.length == 2 and alphabet_counts.values.min > 1)
            newseqs.alignment_concat site
          end #else the column is uninformative and ignored
        end
        
        return newseqs
      end
    end
  end
end



if __FILE__ == $0
  SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')
  
  # Parse command line options into the options hash
  options = {
    :logger => 'stderr',
  }
  o = OptionParser.new do |opts|
    #TODO Fill in usage, description and option parsing below
    opts.banner = "
      Usage: #{SCRIPT_NAME} <aligned_fasta_file>
      
      Remove columns from an alignment that are definitely phylogenetically uninformative. Specifically, remove columns with all gaps, and remove columns where all except 1 character are the same\n\n"
    # Example option
    opts.on("-e", "--eg", "description [default: #{options[:eg]}]") do |f|
      options[:operation] = OVERALL
    end
    
    # logger options
    opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {Bio::Log::CLI.trace('error')}
    opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
    opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| Bio::Log::CLI.trace(s)}
  end
  o.parse!
  if ARGV.length != 1
    $stderr.puts o
    exit 1
  end
  # Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
  Bio::Log::CLI.logger(options[:logger]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)
  
  
  aln = Bio::Alignment.new(Bio::FlatFile.open(ARGV[0]))
  puts aln.trim_uninformative_columns.output_fasta


end #end if running as a script