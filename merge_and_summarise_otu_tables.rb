#!/usr/bin/env ruby
require 'set'
require 'pry'
require 'optparse'
require 'bio-logger'
require 'hopcsv'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
  :summary_level => 5 #family-level
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <arguments>

Merge & collapse OTU tables\n\n"

  opts.on("--list-otu-tables ARG",Array, "Comma-separated list of OTU tables that are just lists of taxonomies without abundance, tab delimited") do |arg|
    options[:list_otu_tables] = arg
  end
  opts.on("--app-otu-tables ARG",Array, "Comma-separated list of APP style OTU tables with taxonomy and abundance, tab delimited") do |arg|
    options[:app_otu_tables] = arg
  end

  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end; o.parse!
if ARGV.length != 0
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)


class Sample
  attr_accessor :otu_abundances, :name
end

# Read in list OTU tables
samples = []
doing_eg = true
options[:list_otu_tables].each do |otu_table|
  log.info "Reading #{otu_table}"
  sample = Sample.new
  sample.name = File.basename otu_table
  sample.otu_abundances = {}
  Hopcsv.foreach(otu_table,"\t") do |row|
    key = row[0...options[:summary_level]]
    if doing_eg
      log.debug "Example lineage: #{key}"
      doing_eg = false
    end
    sample.otu_abundances[key] ||= 0
    sample.otu_abundances[key] += 1
  end
  log.info "Read in #{sample.otu_abundances.length} lineages from #{otu_table}"

  samples.push sample
end

# Read in APP tables
unless options[:app_otu_tables].nil?
  options[:app_otu_tables].each do |otu_table|
    log.info "Reading APP style table #{otu_table}"

    line_number = 0
    name_to_sample = {}
    lineage_column = nil
    doing_eg = true
    Hopcsv.foreach(otu_table,"\t") do |row|
      line_number += 1
      if line_number == 2
        lineage_column = row.find_index 'Consensus Lineage'
        if lineage_column.nil?
          raise "Unable to find lineage column in app OTU table"
        end

        headers = row[1...lineage_column]
        log.info "Found #{headers.length} samples in this OTU table"
        headers.each_with_index do |header|
          name = header
          if name_to_sample.key?(header)
            name = "#{header}:#{(1..100).to_a.sample(1)[0] }"
            log.warn "Renaming #{header} to #{name} to preserve uniqueness"
          end

          name_to_sample[name] = Sample.new
          name_to_sample[name].name = name
          name_to_sample[name].otu_abundances = {}
        end
      end
      next if row[0][0] == '#'

      taxonomy = row[lineage_column].split('; ')[0...options[:summary_level]]
      if doing_eg
        log.debug "Example lineage: #{taxonomy}"
        doing_eg = false
      end
      name_to_sample.keys.each_with_index do |sample, i|
        abundance = row[i+1].to_i
        next if abundance == 0
        name_to_sample[sample].otu_abundances[taxonomy] ||= 0
        name_to_sample[sample].otu_abundances[taxonomy] += abundance
      end
    end

    name_to_sample.each do |name, sample|
      samples.push sample
    end

    log.info "Read #{name_to_sample.length} samples containing e.g. #{name_to_sample[name_to_sample.keys[0]].otu_abundances.length} lineages"
  end
end


# Print out the merged OTU table
lineages = Set.new
samples.each{|s| s.otu_abundances.keys.each {|k| lineages << k}}
log.info "Found #{lineages.length} lineages to report"

# print headers
print "\t"
puts samples.collect{|s| s.name}.join("\t")
# Print data
lineages.each do |lineage|
  print lineage.join('; ')
  print "\t"
  puts samples.collect{|s| s.otu_abundances[lineage].nil? ? 0 : s.otu_abundances[lineage]}.join("\t")
end





