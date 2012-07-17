#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require File.join(File.basename(__FILE__),'454UnscaffoldedContigs')

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')



# Parse command line options into the options hash
options = {
  :logger => 'stderr',
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} 454ContigGraph.txt

    Convert the 454ContigGraph.txt output into a \"scaffoldNumber    scaffoldLength    scaffoldCoverage\" tabular form\n\n"

  # logger options
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {Bio::Log::CLI.trace('error')}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| Bio::Log::CLI.trace(s)}
end
o.parse!
if ARGV.length != 1 and ARGV.length != 0
  $stderr.puts "Unexcepted numebr of arguments found: #{ARGV.length}"
  $stderr.puts o
  exit 1
end
# Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
Bio::Log::CLI.logger(options[:logger]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)

state = 'read_contig_stats'
class Contig
  attr_accessor :length, :coverage
end

contig_stats = {}

puts %w(scaffold length coverage).join("\t")

ARGF.each_line do |line|
  if state == 'read_contig_stats'
    splits = line.strip.split("\t")
    if splits.length == 6
      state = 'scaffolds'
      unless splits[0] == 'C'
        raise "Unexpected line: #{splits.inspect}"
      end
      log.info "Cached stats on #{contig_stats.length} contigs"
    elsif splits.length == 4
      contig = Contig.new
    contig.length = splits[2].to_i
    contig.coverage = splits[3].to_f
    contig_stats[splits[0]] = contig
    else
      raise
    end
  end

  if state == 'scaffolds'
    if line.match(/^S/)
      splits = line.split("\t")
      raise unless splits.length == 4
      total_length = splits[2]

      #caculate coverage, not including gaps into the calculation
      components = splits[3].split(';')
      total_bases = 0
      total_contiged_length = 0
      components.each do |c|
        cees = c.split(':')
        raise unless cees.length == 2
        if cees[0]!='gap'
          contig = contig_stats[cees[0]]
          if contig.nil?
            raise "contig not found: #{cees[0]}"
          end
        total_contiged_length += contig.length
        total_bases += contig.coverage*contig.length
        end
      end
      puts [
        splits[1],
        total_length,
        total_bases/total_contiged_length
      ].join("\t")
    end
  end
end
