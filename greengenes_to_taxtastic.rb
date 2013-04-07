#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'set'
require 'csv'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :output_filenames_base => 'gg',
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} [options] --taxonomy <gg_tax_file>

    Takes a greengenes tax file and outputs a .tax.csv and .seqinfo.csv for use with taxtastic\n\n"

  opts.on("--taxonomy SOME_OTU_TAXONOMY", "taxonomy file from greengenes e.g. 97_otu_taxonomy.txt [required]") do |arg|
    options[:taxonomy_file] = arg
  end
  opts.on("--output-filenames-base BASE", "File name to build upon e.g. 'gg' => './gg.tax.csv' and './gg.seqinfo.csv' [default: #{options[:output_filenames_base]}]") do |arg|
    options[:output_filenames_base] = arg
  end

  # logger options
  opts.separator "\n\tVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {Bio::Log::CLI.trace('error')}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| Bio::Log::CLI.trace(s)}
end; o.parse!
if ARGV.length != 0 or options[:taxonomy_file].nil?
  $stderr.puts o
  exit 1
end
# Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
Bio::Log::CLI.logger(options[:logger]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)

# The tax file must have these 4 fields first:
# tax_id	parent_id	rank	tax_name
tax_output_file = File.open("#{options[:output_filenames_base]}.tax.csv",'w')
seqinfo_output_file = File.open("#{options[:output_filenames_base]}.seqinfo.csv",'w')

# headers
tax_output_file.puts %w(tax_id parent_id rank tax_name).join(',')
seqinfo_output_file.puts %w(seqname tax_id).join(',')

# real data
already_printed_taxids = Set.new
level_to_rank = %w(root
kingdom
phylum
class
order
family
genus
species
)
CSV.foreach(options[:taxonomy_file], :col_sep => "\t") do |row|
  raise unless row.length == 2
  gg_id = row[0]
  tax_string = row[1]
  # 1111750 k__Bacteria; p__TM7; c__TM7-3; o__I025; f__; g__; s__
  taxonomies = tax_string.split('; ')
  raise "Unexpected taxonomies format in line #{row.inspect}" unless taxonomies.length == 7
  raise "Unexpected comma in line #{row.inspect}" if tax_string.match(/,/)
  taxonomies = ['Root']+taxonomies
  # Remove from taxonomies all that are
  first_unclassified_index = taxonomies.index do |taxon|
    taxon.match(/^.__$/)
  end
  unless first_unclassified_index.nil?
    taxonomies = taxonomies[0...first_unclassified_index]
  end

  seqinfo_output_file.puts [
    "gg_otu_#{gg_id}",
    taxonomies.join(';')
  ].join(',')

  breadcrumbs = ''
  taxonomies.each_with_index do |taxon, level|
    new_breadcrumbs = nil
    if level == 0
      new_breadcrumbs = taxon
    else
      new_breadcrumbs = breadcrumbs + ';' + taxon
    end

    unless already_printed_taxids.include?(new_breadcrumbs)
      #tax_output_file.puts [
      puts [
        new_breadcrumbs,
        level == 0 ? 'Root' : breadcrumbs, #Root is special, its parent is itself
        level_to_rank[level],
        taxon.gsub(/^.__/,''),
      ].join(',')
      already_printed_taxids << new_breadcrumbs
    end
    breadcrumbs = new_breadcrumbs
  end
end



tax_output_file.close
seqinfo_output_file.close













