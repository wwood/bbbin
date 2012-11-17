#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'csv'

if __FILE__ == $0 #needs to be removed if this script is distributed as part of a rubygem
  SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')
  
  # Parse command line options into the options hash
  options = {
    :logger => 'stderr',
  }
  o = OptionParser.new do |opts|
    opts.banner = "
      Usage: #{SCRIPT_NAME} <arguments>
      
      Take a tab separated file of sample IDs and more useful names, together with an APP (ACE Pytotag Pipeline) config file.
      Replace the sample IDs in the config file with the ones from the tab-separated file, for easier reading. \n\n"
      
    opts.on("-c", "--config APP_CONFIG_FILE", "ACE Pyrotag Pipeline config file [required]") do |arg|
      options[:config_file_path] = arg
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
  if ARGV.length != 0 or options[:config_file_path].nil? or options[:rename_key_path].nil?
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
  
  # Read the config file and replace each of the things
  state = 'top'
  CSV.foreach(options[:config_file_path],:col_sep => "\t") do |row|
    if state=='top'
      if row[0]=='@@'
        puts row.join("\t")
        state = 'bottom'
      elsif row[0].strip[0]=='#'
        puts row.join("\t")
      else
        current = row[0]
        if renames.key?(current)
          row[0] = renames[current]
        else
          log.warn "No rename key found for #{current}, so just letting this sample pass through untouched"
        end
        puts row.join("\t")
      end
    elsif state=='bottom'
      puts row.join("\t")
    end
  end
end #end if running as a script