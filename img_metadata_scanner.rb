#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'

$:.unshift File.join(ENV['HOME'],'git','bioruby-taxonomy_definition_files','lib')
require 'bio-taxonomy_definition_files'

if __FILE__ == $0 #needs to be removed if this script is distributed as part of a rubygem
  SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')
  IMG_METADATA_FILE_ENV_VARIABLE_NAME='IMG_METADATA_FILE'
  
  # Parse command line options into the options hash
  options = {
    :logger => 'stderr',
    :output_fields => ['taxon_oid'],
    :listing_mode => false,
  }
  o = OptionParser.new do |opts|
    opts.banner = "
      Usage: #{SCRIPT_NAME} -i <img_metadata_file> <key=value>
      
      Runs through an IMG metadata file, and prints out the IMG identifiers of those genomes that match the key=value criteria. E.g. Kingdom=Archaea to grep for all archaeons.\n\n"
      
    opts.on("-i", "--img-metadata-file PATH", "Path to IMG metadata file [required]. This is not necessary if there is a valid environment variable #{IMG_METADATA_FILE_ENV_VARIABLE_NAME} available.") do |arg|
      options[:img_metadata_file] = arg
    end
    opts.on("-o", "--output-fields FIELDS", "List of output fields, comma separated [default: #{options[:output_fields].join(',')}]") do |arg|
      options[:output_fields] = arg.split ','
    end
    opts.on("-l", "--list", "Instead of filtering print a list of the fields in the output file, newline separated [default: #{options[:listing_mode]}]") do |arg|
      options[:listing_mode] = true
    end

    # logger options
    opts.separator "\nVerbosity:\n\n"
    opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {Bio::Log::CLI.trace('error')}
    opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
    opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| Bio::Log::CLI.trace(s)}
  end; o.parse!

  # Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
  Bio::Log::CLI.logger(options[:logger]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)
  
  # Read in the system-specific ENV variable for the path to the img metadata file unless a path has already been specified 
  if options[:img_metadata_file].nil?
    if ENV[IMG_METADATA_FILE_ENV_VARIABLE_NAME].nil?
      $stderr.puts "img-metadata-file not found as either a command line option or environment variable, failing."
      $stderr.puts o
      exit 1
    else
      options[:img_metadata_file] = ENV[IMG_METADATA_FILE_ENV_VARIABLE_NAME]
      log.debug "Using environment variable #{IMG_METADATA_FILE_ENV_VARIABLE_NAME} to define path to IMG metadata file #{options[:img_metadata_file]}"
    end
  end
  unless File.exist?(options[:img_metadata_file])
    $stderr.puts "IMG metadata file #{options[:img_metadata_file]} not found - was it specified correctly?"
    exit 2
  end
  

  
  # If listing mode is set, don't do any filtering, just list the variables
  if options[:listing_mode]
    # Read the IMG metadata
    img = Bio::IMG::TaxonomyDefinitionFile.read(options[:img_metadata_file])
    log.info "Found #{img.length} taxons in the IMG metadata file"
    
    img[0].attributes.keys.each do |key|
      puts key
    end
    
    exit 0
  end
    
  # Parse the key/value pairs
  filter_hash = {}
  unless ARGV[0].nil?
    ARGV[0].split(',').each do |split|
      splits2 = split.split('=')
      unless splits2.length == 2
        log.error "Badly parsed key/value pair '#{split}', expected exactly 1 = sign, found #{splits2.length}"
        exit 1
      end
      
      key = splits2[0]
      unless filter_hash[key].nil?
        log.error "Duplicate filter key found: #{key}, failing"
        exit 1
      end
      filter_hash[key] = splits2[1]
    end
  end
  log.info "Using #{filter_hash.length} filters."
  
  # Read the IMG metadata
  img = Bio::IMG::TaxonomyDefinitionFile.read(options[:img_metadata_file])
  log.info "Found #{img.length} taxons in the IMG metadata file"
  
  # Expect that each of the filters are found in the list of headers available, otherwise filtering will be ineffectual
  # Ditto for the output names
  [filter_hash.keys, options[:output_fields]].flatten.each do |key|
    unless img[0].attributes.keys.include?(key)
      log.warn "Unable to find column named #{key} in the IMG metadata file - typo perhaps?"
      exit 1
    end
  end
  
  # Go through each row, printing the outputs if they pass the filter
  img.each do |taxon|
    passed = true
    filter_hash.each do |key, value|
      if taxon.attributes[key] != value
        passed = false
        break
      end
    end
    
    if passed
      puts options[:output_fields].collect{|field| taxon.attributes[field]}.join("\t")
    end
  end
  
end #end if running as a script