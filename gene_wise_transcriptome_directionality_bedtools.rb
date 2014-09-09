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


calculate_cov = lambda do |covs|
  covs.reduce(:+).to_f/covs.length
end

get_covs = lambda do |cov_lines|
  feature_to_covs = {}
  previous_feature = nil
  covs = []
  cov_lines.each_line do |line|
    splits = line.split("\t")
    #gi|169887498|gb|CP000948.1|
    #Prodigal_v2.6.1
    #CDS
    #1047994
    #1049139
    #157.2
    #-
    #0
    #ID=1_972;partial=00;start_type=ATG;rbs_motif=AGGA;rbs_spacer=5-10bp;gc_cont=0.568;conf=100.00;score=157.25;cscore=141.04;sscore=16.20;rscore=10.98;uscore=-0.89;tscore=3.93;
    #1047994
    #111
    feat = splits[0..8]
    if feat != previous_feature
      feature_to_covs[previous_feature] = calculate_cov.call(covs) unless previous_feature.nil?
      covs = []
    end
    covs.push splits[10].to_i
    previous_feature = feat
  end
  feature_to_covs[previous_feature] = calculate_cov.call(covs)

  feature_to_covs
end



cov_lines_fwd = Bio::Commandeer.run "bedtools coverage -abam #{bam_file.inspect} -b #{gff_file.inspect} -d -s", :log => log
#cov_lines_fwd = Bio::Commandeer.run "samtools view -b #{bam_file.inspect} 'gi|169887498|gb|CP000948.1|:1-2000' |bedtools coverage -abam - -b #{gff_file.inspect} -d -s", :log => log
log.info "Parsing forward aligned reads"
covs_fwd = get_covs.call(cov_lines_fwd)

cov_lines_rev = Bio::Commandeer.run "bedtools coverage -abam #{bam_file.inspect} -b #{gff_file.inspect} -d -S", :log => log
#cov_lines_rev = Bio::Commandeer.run "samtools view -b #{bam_file.inspect} 'gi|169887498|gb|CP000948.1|:1-2000' |bedtools coverage -abam - -b #{gff_file.inspect} -d -S", :log => log
log.info "Parsing reverse aligned reads"
covs_rev = get_covs.call(cov_lines_rev)

puts [
  'contig',
  'type',
  'start',
  'end',
  'strand',
  'forward_average_coverage',
  'reverse_average_coverage',
  'annotation',
].join("\t")

covs_fwd.each do |feature, cov_fwd|
  cov_rev = covs_rev[feature]
  record = Bio::GFF::GFF3::Record.new(feature.join("\t"))

  products = record.attributes.select{|a| a[0] == 'product'}
  product = 'unannotated'
  if products.length == 1
    product = products[0][1]
  end

  puts [
    record.seqname,
    record.feature,
    record.start,
    record.end,
    record.strand,
    cov_fwd,
    cov_rev,
    product,
  ].join("\t")
end
