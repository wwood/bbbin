#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'systemu'
require 'tempfile'
require 'bio'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} -c contigs.list -f fasta.fa

    Given a list of contigs and a fasta file containing all of those contigs (+possibly others), return a bed file of those contigs encompassing their whole length\n\n"

  opts.on("-c", "--contig-list FILE", "list of contigs, newline separated [required]") do |arg|
    options[:contig_list] = arg
  end
  opts.on("-f", "--fasta FILE", "Fasta file of contigs [required]") do |arg|
    options[:fasta_file] = arg
  end

  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end; o.parse!
if ARGV.length != 0 or options[:contig_list].nil? or options[:fasta_file].nil?
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)

contig_list = []
File.open(options[:contig_list]).each_line do |line|
  name = line.strip
  #log.debug "Adding #{name} to contig list"
  contig_list.push name
end

list_string = contig_list.collect{|c| c.inspect}.join(' ')
Tempfile.open(SCRIPT_NAME) do |tempfile|
  faidx_command = "samtools faidx #{options[:fasta_file].inspect} #{list_string} >#{tempfile.path}"
  log.debug "Running command: #{faidx_command}"
  status, stdout, stderr = systemu faidx_command
  raise "Error running samtools: #{stderr}" unless status.exitstatus == 0


  num_recovered = 0
  Bio::FlatFile.foreach(tempfile.path) do |seq|
    num_recovered +=1
    puts [
      seq.definition,
      1,
      seq.seq.length
    ].join("\t")
  end
  if num_recovered == contig_list.length
    log.info "Successfully created bed file for #{num_recovered} contigs"
  else
    log.error "#{contig_list.length-num_recovered} contigs not converted to bed file properly! Do they exist in the fasta file?"
    exit 1
  end
end

