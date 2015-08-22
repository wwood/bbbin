#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <arguments>

    Description of what this program does...\n\n"

  opts.on("--taxonomy", "IMG taxonomy, tab separated ID then taxonomy [required]") do |arg|
    options[:taxonomy] = arg
  end
  opts.on("--type-strains", "type strain list, tab separated ID then taxonomy [required]") do |arg|
    options[:type_strains] = arg
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
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME); log.outputters[0].formatter = Log4r::PatternFormatter.new(:pattern => "%5l %c %d: %m", :date_pattern => '%d/%m %T')


# read type strains, hash by taxonomy
type_strain_hash = {}
CSV.foreach(options[:type_strains], :col_sep => "\t") do |row|
  raise unless row.length == 2
  type_strain_hash[row[1]] ||= []
  type_strain_hash[row[1]].push row[2]
end
log.info "Read in #{type_strain_hash.length} different species of type strains"
# read taxonomy, hash by taxonomy
taxonomy_hash = {}
CSV.foreach(options[:taxonomy], :col_sep => "\t") do |row|
  raise unless row.length == 2
  taxonomy_hash[row[1]] ||= []
  taxonomy_hash[row[1]].push row[2]
end
log.info "Read in #{taxonomy_hash.length} different species in total"

# for each tax
taxonomy_hash.each do |tax, img_ids|
  # accept if it is s__
  if img_ids.match(/s__$/)
    img_ids.each do |img|
      puts [img, tax].join("\t")
    end
  elsif type_strain_hash.key?(tax)# choose a random type strain if there is one
    puts [type_strain_hash[tax].sample, tax].join("\t")
  else# else choose a random one from the list
    puts [img_ids.sample, tax].join("\t")
  end
end
