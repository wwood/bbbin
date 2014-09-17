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
  opts.on("--require-start", "if set, require mapped reads to start within the gene boundary (not just overlap with the gene) in order for them to count towards the gene's total [default: don't require]") do
    options[:require_start] = true
  end
  opts.on("--single", "assume reads mapped were single ended [default: assume paired-end]]") do
    options[:single_ended] = true
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

headers = [
  'contig',
  'type',
  'start',
  'end',
  'strand',
]
fwd_flag = 'pPr1'
rev_flag = 'pPR1'
if options[:single_ended]
  fwd_flag = ''
  rev_flag = 'r'
  headers.push 'num_forward'
  headers.push 'num_reverse'
else
  headers.push 'num forward read 1'
  headers.push 'num reverse read 1'
end
all_flags = [fwd_flag, rev_flag]
headers.push 'FPKG'
headers.push 'annotation'
puts headers.flatten.join("\t")

num_genes_printed = 0
view_params = '-f 66 -F3336'
view_params = '-F 3841' if options[:single_ended]

Bio::FlatFile.open(Bio::GFF::GFF3, gff_file).entries[0].records.each do |record|
  next if record.start.nil? #Ignore comment lines etc.

  position = "#{record.seqname}:#{record.start}-#{record.end}"
  sams = Bio::Commandeer.run "samtools view -X #{view_params} #{bam_file.inspect} #{position.inspect}"
  flags_hash = {}
  seen_reads = Set.new
  gene_start = record.start
  sams.each_line do |sam|
    r = sam.split("\t")
    read = r[0]
    flags = r[1]
    start = r[3].to_i

    if !options[:require_start] or start >= gene_start
      flags_hash[flags] ||= 0
      flags_hash[flags] += 1
    end
  end

  products = record.attributes.select{|a| a[0] == 'product'}
  product = 'unannotated'
  if products.length == 1
    product = products[0][1]
  end

  num_fwd = flags_hash[fwd_flag]; num_fwd ||= 0
  num_rev = flags_hash[rev_flag]; num_rev ||= 0
  fpkg = (num_fwd+num_rev).to_f/(record.end-record.start)*1000

  puts [
    record.seqname,
    record.feature,
    record.start,
    record.end,
    record.strand,
    num_fwd,
    num_rev,
    fpkg,
    product,
  ].join("\t")

  num_genes_printed += 1
  break if options[:max_genes] and num_genes_printed > options[:max_genes]
end
