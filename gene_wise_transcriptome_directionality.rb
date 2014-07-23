#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'bio'
require 'bio-commandeer'

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

  opts.on("-e", "--eg ARG", "description [default: #{options[:eg]}]") do |arg|
    options[:example] = arg
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

gff_file = ARGV[0]
bam_file = ARGV[1]


parsed_records = []
Bio::FlatFile.open(Bio::GFF::GFF3, ARGV[0]).entries[0].records.each do |record|
  sams = Bio::Commandeer.run "samtools view #{bam_file.inspect} #{record.seqname}#{record.start}-#{record.end}"
flags_hash = {}
sams.each do |sam|
  r = sam.split("\t")
  flags = r[1]
  sams[flags] ||= 0
  sams[flags] += 1
end

parsed_records.push [
  flags_hash,
  record.seqname,
  record.feature,
  record.start,
  record.end,
record.strand,
  ]
end

all_flags = Set.new
parsed_records.each do |r|
  r[0].keys.each do |f|
    all_flags << f
  end
end
flags = all_flags.to_a.sort

# headers
puts [
  flags,
  'contig',
  'type',
  'start',
  'end',
  'strand'
  ].flatten.join("\t")

parsed_records.each do |r|
  print r[1..-1].join("\t")
  print "\t"
  all_flags.each do |f|
    print "\t"
    if r[0].key?(f)
      print r[0][f]
    else
      print '0'
    end
  end
  puts
end
