#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'bio-cigar'
require 'bio-samtools'
require 'bio'
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

    Takes a SAM file on stdin and prints a line for each read: hit_id, %ID, num_matching, num_not_matching (tab-separated).

    Intended to be used with a samfile without headers piped in (i.e. the default output from 'samtools view -F4 <bam>')
    \n\n"

  opts.on("-r", "--reference-seqs PATH", "path to a fasta file of the reference sequences refered to in the samfile [required]") do |arg|
    options[:seqfile] = arg
  end
  opts.on("--first-hit-only", "Only print the first hit of each query (sometimes sam files contain multiple lines for the same read) [default: no]") do
    options[:first_hit_only] = true
  end
  opts.on("--min-length LENGTH", "Don't print when the match length is less than this length [default: no length restriction]") do |arg|
    options[:min_length] = arg.to_i
  end
  opts.on("--min-identity", "--min-percent PERCENT", "Don't print when the percent identity of the match is less than this length [default: no restriction]") do |arg|
    options[:min_percent] = arg.to_i
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

# Read in fasta seqs
seqs = {}
log.info "Reading reference sequence file.."
Bio::FlatFile.foreach(options[:seqfile]) do |seq|
  seqs[seq.definition.split(/\s+/)[0]] = seq.seq
end
log.info "Read in #{seqs.length} sequences from the reference sequence file."

count = 0
already_seen_qnames = Set.new
ARGF.each_line do |line|
  aln = Bio::DB::Alignment.new
  aln.sam = line

  if options[:first_hit_only]
    next if already_seen_qnames.include?(aln.qname)
    already_seen_qnames << aln.qname
  end

  ref = seqs[aln.rname]
  raise "Couldn't find sequence #{aln.rname}" if ref.nil?
  begin
    percent, pos, neg = aln.percent_identity(ref)
    next if options[:min_length] and pos+neg<options[:min_length]
    next if options[:percent_identity] and percent < options[:percent_identity]
    puts [
      aln.qname,
      aln.rname,
      percent,
      pos,
      neg,
      aln.seq
    ].join("\t")
    count += 1
  rescue Exception => e
    log.warn e.to_s
  end
end
log.info "Printed information on #{count} alignments"
