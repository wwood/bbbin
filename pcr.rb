#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'bio-ipcress'
require 'bio'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
  :print_amplicon => false,
  :num_mismatches => 3,
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
  opts.on("--fasta FASTA_FILE[, FASTA_FILE2, ..]", Array, "sequence(s) being assayed [required]") do |arg|
    options[:fasta] = arg
  end
  opts.on("--mismatches NUM", Integer, "max number of allowed mismatches. [default: #{options[:num_mismatches] }]") do |arg|
    options[:num_mismatches] = arg
  end
  opts.on("--print-amplicon", "print out the sequence of the amplicon [default: #{options[:print_amplicon] }]") do
    options[:print_amplicon] = true
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
headers = %w(
target
mismatches_fwd
mismatches_rev
length
gc_of_forward_matching
gc_of_reverse_matching
gc_positions_of_forward_matching
gc_positions_of_reverse_matching
)
headers += ['amplicon'] if options[:print_amplicon]
puts headers.join("\t")

options[:fasta].each do |fasta|
  # run ipcress
  mismatch_param = 3 #default to 3 so ipcress bug gets worked around, and filter later
  mismatch_param = options[:num_mismatches] if options[:num_mismatches] > mismatch_param

  results = Bio::Ipcress.run primer_set, fasta, :mismatches => mismatch_param

  seqs = {}
  if options[:print_amplicon]
    Bio::FlatFile.foreach(fasta) do |e|
      name = e.definition
      seq = e.seq.seq
      raise "Duplicate sequence name found: #{name}" if seqs.key?(name)
      seqs[name] = seq
    end
  end

  results.each do |res|
    misses = res.recalculate_mismatches_from_alignments
    next if misses.reduce(:+) > options[:num_mismatches]

    to_print = [
      res.target,
      misses[0],
      misses[1],
      res.length,
      to_gc_count.call(res.forward_matching_sequence),
      to_gc_count.call(res.reverse_matching_sequence),
      to_gc_binary.call(res.forward_matching_sequence),
      to_gc_binary.call(res.reverse_matching_sequence),
      ]
    if options[:print_amplicon]
      name = res.target.gsub(':filter(unmasked)','') #H509DRAFT_scaffold00021.21:filter(unmasked)
      raise "Unable to find sequence name #{name} in fasta file, possible programming error" if !seqs.key?(name)
      seq = seqs[name]
      amplicon = seq[(res.start)...(res.start+res.length)]
      if res.result_type == 'forward'
        to_print += [amplicon]
      elsif res.result_type == 'revcomp'
        to_print += [Bio::Sequence::NA.new(amplicon).reverse_complement.to_s.upcase]
      else
        raise "Unexpected ipcress result type: #{res.result_type}"
      end
    end
    puts to_print.join("\t")
  end
end
