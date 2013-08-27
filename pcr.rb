#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'bio-ipcress'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <arguments>

    Take a set of primers and a fasta file of genome(s), and output each place that the primer hits \n\n"

  opts.on("--primer1 PRIMER", "sequence of the forward primer [required]") do |arg|
    options[:primer1] = arg
  end
  opts.on("--primer2 PRIMER", "sequence of the reverse primer [required]") do |arg|
    options[:primer2] = arg
  end
  opts.on("--fasta FASTA_FILE", "sequence of the genome(s) being assayed [required]") do |arg|
    options[:fasta] = arg
  end

  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end; o.parse!
if ARGV.length != 0 or options[:primer1].nil? or options[:primer2].nil? or options[:fasta].nil?
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)


# make the primer set
primer_set = Bio::Ipcress::PrimerSet.new options[:primer1], options[:primer2]

# run ipcress
results = Bio::Ipcress.run primer_set, options[:fasta], :mismatches => 3

to_gc_binary = lambda do |seq|
  str = ''
  seq.each_char do |char|
    if %(G C).include?(char)
      str="#{str}1"
    else
      str="#{str}0"
    end
  end
  str
end
to_gc_count = lambda do |seq|
  count = 0
  seq.each_char do |char|
    count += 1 if %(G C).include?(char)
  end
  count
end

# output characters of each hit
puts %w(
target
mismatches_fwd
mismatches_rev
length
gc_of_forward_matching
gc_of_reverse_matching
gc_positions_of_forward_matching
gc_positions_of_reverse_matching
).join("\t")
results.each do |res|
  misses = res.recalculate_mismatches_from_alignments

  puts [
    res.target,
    misses[0],
    misses[1],
    res.length,
    to_gc_count.call(res.forward_matching_sequence),
    to_gc_count.call(res.reverse_matching_sequence),
    to_gc_binary.call(res.forward_matching_sequence),
    to_gc_binary.call(res.reverse_matching_sequence),
  ].join("\t")
end
