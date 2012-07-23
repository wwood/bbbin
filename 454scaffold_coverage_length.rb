#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require File.join(File.dirname(__FILE__),'454UnscaffoldedContigs')

$:.unshift File.join([ENV['HOME'], %w(git bioruby-agp lib)].flatten)
require 'bio-agp'

$:.unshift File.join([ENV['HOME'], %w(git bioruby-newbler_outputs lib)].flatten)
require 'bio-newbler_outputs'

script_name = File.basename(__FILE__); log_name = script_name.gsub('.rb','')



# Parse command line options into the options hash
options = {
  :logger => 'stderr',
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{script_name} <arguments>

    Convert the newbler outputs into a \"scaffoldNumber    scaffoldLength    scaffoldCoverage\" tabular form\n\n"

  opts.on("-a", "--alignment-info 454ALIGNMENT_INFO.TXT", "From newbler [required unless -c/--contig-graph is defined]") do |f|
    options[:alignment_info_file] = f
  end
  opts.on("-s", "--scaffolds-txt 454SCAFFOLDS.TXT", "From newbler [required unless -c/--contig-graph is defined]") do |f|
    options[:scaffolds_file] = f
  end
  opts.on("-c", "--contig-graph CONTIG_GRAPH.TXT", "From newbler [required unless both -a and -s are defined]") do |f|
    options[:contigs_graph_file] = f
  end

  # logger options
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {Bio::Log::CLI.trace('error')}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| Bio::Log::CLI.trace(s)}
end
o.parse!
if ARGV.length != 1 and ARGV.length != 0
  $stderr.puts "Unexpected number of arguments found: #{ARGV.length}"
  $stderr.puts o
  exit 1
end
# Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
Bio::Log::CLI.logger(options[:logger]); log = Bio::Log::LoggerPlus.new(log_name); Bio::Log::CLI.configure(log_name)

state = 'read_contig_stats'
class Contig
  attr_accessor :length, :coverage
end

contig_stats = {}

puts %w(scaffold length coverage).join("\t")


# If using contig_graph.txt
if !options[:contigs_graph_file].nil?
  File.open(options[:contigs_graph_file]).each_line do |line|
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
  
  
# Using alignment file and scaffolds (AGP) file
elsif options[:alignment_info_file] and options[:scaffolds_file]
  log.info "Caching contig lengths and coverages from #{options[:alignment_info_file]}"
  contig_hash = Bio::Newbler::AlignmentInfoFile.contig_hash(options[:alignment_info_file])
  log.info "Cached contig information for #{contig_hash.length} contigs"
  #require 'pp'
  #pp contig_hash
  
  log.info "Iterating AGP file #{options[:scaffolds_file]}"
  Bio::Assembly::AGP.new(options[:scaffolds_file]).each_scaffold do |scaffold|
    contig_coverages = []
    contig_lengths = []
    total_scaffold_length = 0
    
    scaffold.components.each_with_index do |component, i|
      # Add the length of the component
      total_scaffold_length += component.length
      
      # Don't count gaps in the coverages
      next unless component.kind_of?(Bio::Assembly::AGP::Scaffold::Object)
      
      contig_name = component.component_id
      contig = contig_hash[contig_name]
      #contig = contig_hash[contig_hash.keys[i]]
      
      if contig.nil?
        raise "Unable to find contig #{contig_name} in alignment info file #{options[:alignment_info_file]}"
      end
      
      contig_coverages.push contig.coverage
      contig_lengths.push component.length #maybe not all the contig is included in the scaffold? Use component length to be safe
    end
    
    #coverage (coverage1*length1+coverage2*length2)/(length1+length2)
    numerator = 0
    contig_lengths.each_with_index do |len, i|
      numerator += contig_coverages[i]*len
    end
    coverage = numerator/contig_lengths.inject{|sum,i|sum+=i}
    
    # Print out the answers
    puts [
      scaffold.identifier,
      total_scaffold_length,
      coverage,
    ].join("\t")
  end
  


# Else total fail
else
  $stderr.puts o
  exit 1
end
