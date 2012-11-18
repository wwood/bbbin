#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'csv'

$:.push File.join(File.dirname(__FILE__),'..','bioruby-otu_table','lib')
require 'bio-otu_table'

if __FILE__ == $0 #needs to be removed if this script is distributed as part of a rubygem
  SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = 'bio-otu_table'
  
  # Parse command line options into the options hash
  options = {
    :logger => 'stderr',
  }
  o = OptionParser.new do |opts|
    opts.banner = "
      Usage: #{SCRIPT_NAME} <arguments>
      
      Take a tab separated file of sample IDs and more useful names, together with an OTU table.
      Replace the sample IDs in the OTU table with the ones from the tab-separated file, for easier reading. \n\n"
      
    opts.on("-i", "--otu-table OTU_TABLE_FILE", "OTU table to be manipulated [required]") do |arg|
      options[:otu_table_path] = arg
    end
    opts.on("-t", "-k", "--rename-key TAB_SEPARATED_RENAME_KEY_FILE", "File of old names and new names (in that order, no headers [required]") do |arg|
      options[:rename_key_path] = arg
    end

    # logger options
    opts.separator "\nVerbosity:\n\n"
    opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {Bio::Log::CLI.trace('error')}
    opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
    opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| Bio::Log::CLI.trace(s)}
  end; o.parse!
  if ARGV.length != 0 or options[:otu_table_path].nil? or options[:rename_key_path].nil?
    $stderr.puts o
    exit 1
  end
  # Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
  Bio::Log::CLI.logger(options[:logger]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)
  
  
  # Read keys file
  renames = {}
  CSV.foreach(options[:rename_key_path],:col_sep => "\t") do |row|
    old = row[0]
    newer = row[1]
    unless renames[old].nil?
      raise "Duplicate key found in ID file, fail: #{old}"
    end
    renames[old] = newer
  end
  
  # Read the OTU table
  otu_table = Bio::OtuTable.read(options[:otu_table_path])
  
  # Rename, keeping the order. In Ruby, hashes remain ordered.
  old_keys = otu_table.sample_names
  new_otu_table = Bio::OtuTable.new
  new_otu_table.samples = {}
  new_otu_table.otu_identifiers = otu_table.otu_identifiers
  
  log.debug "Before renaming, there is #{otu_table.samples.length} samples"
  otu_table.sample_names.each do |sample_name|
    if renames.key? sample_name
      new_name = renames[sample_name]
      values = otu_table.samples[sample_name]
      log.debug "Renaming sample name #{sample_name} to #{new_name}"
      new_otu_table.samples[new_name] = values
    else
      log.warn "Unable to rename sample called #{sample_name} since this was not specified in the rename key file"
      new_otu_table.samples[sample_name] = otu_table.samples[sample_name]
    end
  end
  log.debug "After renaming, there is still #{new_otu_table.samples.length} samples"
  print new_otu_table.to_s
end #end if running as a script