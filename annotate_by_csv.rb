#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'csv'
require 'bio'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <arguments>

    Take a blast output and a csv file containing annotations, and a fasta file of the query sequences, and annotate the fasta sequence with the csv annotations.\n\n"

  opts.separator "\nRequired:\n\n"
  opts.on("--annotation FILE", "tab separated annotation file (one column hit id, second is annotation)[required]") do |arg|
    options[:annotation_file] = arg
  end
  opts.on("--blast-result FILE", "outfmt 6 blast output, where the hit IDs correspond to the annotation tsv key column [required]") do |arg|
    options[:blast_file] = arg
  end
  opts.on("--query-fasta FILE", "tab separated annotation file [required]") do |arg|
    options[:query_fasta_file] = arg
  end

  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end; o.parse!
if ARGV.length != 0 or options[:annotation_file].nil? or options[:blast_file].nil? or options[:query_fasta_file].nil?
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)


# Read annotation CSV into hash
log.info "Reading annotations"
hit_to_annotation = {}
eg_hit_id = nil
CSV.foreach(options[:annotation_file], :col_sep => "\t") do |row|
  raise "Unexpected line found in annotation file, I expect ID 'tab annotation' on each line. Found:\n#{row.inspect}" if row.length != 2

  hit_id = row[0]
  annotation = row[1]

  eg_hit_id ||= hit_id

  if hit_to_annotation.key?(hit_id)
    log.warn "Duplicate annotation found, ignoring all but the first, for #{hit_id}"
    next
  end
  hit_to_annotation[hit_id] = annotation
end
log.info "Read in #{hit_to_annotation.length} annotations e.g. #{eg_hit_id} => #{hit_to_annotation[eg_hit_id]}"


# Read blast results into hash
log.info "Reading BLAST results"
query_to_hit = {}
eg_query_id = nil
CSV.foreach(options[:blast_file], :col_sep => "\t") do |row|
  raise "Unexpected line found in blast file, wrong number of columns (#{row.length})" if row.length != 12

  query_id = row[0]
  hit_id = row[1]

  eg_query_id ||= query_id

  if query_to_hit.key?(query_id)
    log.warn "Duplicate blast hit found, ignoring all but the first, for #{query_id}"
    next
  end
  query_to_hit[query_id] = hit_id
end
log.info "Read in #{query_to_hit.length} annotations e.g. #{eg_query_id} => #{query_to_hit[eg_query_id]}"



# Read and annotate the fasta file on the fly
log.info "Writing out annotations of the query sequence file.."
# Annotate the fasta file
annotated_count = 0
no_blast_hit_count=0
no_annotation_count=0
Bio::FlatFile.foreach(options[:query_fasta_file]) do |seq|
  firstname = seq.definition.gsub(/\s.*/,'')
  if query_to_hit.key?(firstname)
    annotated_count += 1
    hit_id = query_to_hit[firstname]
    annotation = hit_to_annotation[hit_id]
    if annotation
      puts ">#{firstname} #{annotation}"
    else
      no_annotation_count += 1
      log.warn "Unable to find annotation for #{hit_id}" unless hit_to_annotation[hit_id]
      puts ">#{firstname}"
    end
  else
    no_blast_hit_count += 1
    puts ">#{firstname}"
  end
  puts seq.seq
end
log.info "Output #{annotated_count} annotated genes, and #{no_blast_hit_count+no_annotation_count} unannotated ones. #{no_blast_hit_count} had no blast and #{no_annotation_count} ones that had no annotation"
