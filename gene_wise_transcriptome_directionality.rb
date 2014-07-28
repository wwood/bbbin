#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'bio'
require 'bio-commandeer'
require 'pry'
require 'set'

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

  opts.on("--bam FILE", "path to mapping file [required]") do |arg|
    options[:bam] = arg
  end
  opts.on("--gff FILE", "path to GFF3 file [required]") do |arg|
    options[:gff] = arg
  end
  opts.separator "\nOptions:\n\n"
  opts.on("--max-genes NUM", Integer, "only parse this many genes [default: not used]") do |arg|
    options[:max_genes] = arg
  end

  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end; o.parse!
if ARGV.length != 0 or options[:bam].nil? or options[:gff].nil?
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)

gff_file = options[:gff]
bam_file = options[:bam]


parsed_records = []

Bio::FlatFile.open(Bio::GFF::GFF3, gff_file).entries[0].records.each do |record|
  sams = Bio::Commandeer.run "samtools view -X -f 66 -F3336 #{bam_file.inspect} #{record.seqname}:#{record.start}-#{record.end}"
  flags_hash = {}
  seen_reads = Set.new
  sams.each_line do |sam|
    r = sam.split("\t")
    read = r[0]
    flags = r[1]

    flags_hash[flags] ||= 0
    flags_hash[flags] += 1
  end

  parsed_records.push [
    flags_hash,
    record.seqname,
    record.feature,
    record.start,
    record.end,
    record.strand,
  ]
  break if options[:max_genes] and parsed_records.length > options[:max_genes]
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
  'contig',
  'type',
  'start',
  'end',
  'strand',
  flags,
  ].flatten.join("\t")

parsed_records.each do |r|
  print r[1..-1].join("\t")
  flags.each do |f|
    print "\t"
    if r[0].key?(f)
      print r[0][f]
    else
      print '0'
    end
  end
  puts
end
