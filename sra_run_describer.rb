#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'

require 'bio-sra'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} -d <SRAmetadb.sqlite> -r <run_list>

    Print descriptions of SRA run accessions\n\n"

  opts.on("-d", "--db PATH", "Path to the SRAmetadb.sqlite file/database [required]") do |f|
    options[:db_path] = f
  end
  opts.on("-r", "--run-file PATH", "Path to a list of SRA run accessions [required]") do |f|
    options[:run_list_path] = f
  end

  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end; o.parse!
if ARGV.length != 0 or options[:db_path].nil? or options[:run_list_path].nil?
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)

include Bio::SRA::Tables

# Connect
Bio::SRA::Connection.connect(options[:db_path])
log.info "Connected to SRAdb database #{options[:db_path]}"

num_described = 0
puts SRA.content_columns.collect{ |cc| cc.name }.join("\t")
File.foreach(options[:run_list_path]) do |line|
  run_id = line.strip
  s = SRA.where(:run_accession => run_id).first
  if s
    puts SRA.content_columns.collect{ |cc|
      r = s.send(cc.name.to_sym)
      r = r.to_s.gsub("\n"," ") unless r.nil?
      r
    }.join("\t")
    num_described += 1
  else
    log.warn "Unable to find SRA accession '#{run_id}'"
  end
end

log.info "Printed descriptions for #{num_described} runs."
