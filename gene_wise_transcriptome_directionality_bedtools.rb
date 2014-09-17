#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'bio'
require 'bio-commandeer'
#require 'pry'
require 'set'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :ignore_directions => false,
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
  opts.separator "\nOptional parameters:\n\n"
  opts.on("--ignore-directions", "ignore directionality, give overall coverage [default: false i.e. differentiate between directions]") do |arg|
    options[:ignore_directions] = true
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
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME); log.outputters[0].formatter = Log4r::PatternFormatter.new(:pattern => "%5l %c %d: %m", :date_pattern => '%d/%m %T')

gff_file = options[:gff]
bam_file = options[:bam]


calculate_cov = lambda do |covs, num_covs|
  covs.reduce(:+).to_f/num_covs
end

get_covs = lambda do |cov_lines|
  feature_to_covs = {}
  previous_feature = nil
  covs = []
  num_covs = 0
  cov_lines.each_line do |line|
    splits = line.split("\t")
    break if splits[0] == 'all'

    #gi|169887498|gb|CP000948.1|
    #Prodigal_v2.6.1
    #CDS
    #1047994
    #1049139
    #157.2
    #-
    #0
    #ID=1_972;partial=00;start_type=ATG;rbs_motif=AGGA;rbs_spacer=5-10bp;gc_cont=0.568;conf=100.00;score=157.25;cscore=141.04;sscore=16.20;rscore=10.98;uscore=-0.89;tscore=3.93;
    #70 #coverage
    #2 #num reads with coverage 70
    #96 #coverage
    #0.0208333
    feat = splits[0..8]
    if feat != previous_feature
      feature_to_covs[previous_feature] = calculate_cov.call(covs, num_covs) unless previous_feature.nil?
      covs = []
      num_covs = 0
    end
    num = splits[10].to_i
    covs.push num*splits[9].to_i
    num_covs += num
    previous_feature = feat
  end
  feature_to_covs[previous_feature] = calculate_cov.call(covs, num_covs)

  feature_to_covs
end


cmd1 = "bedtools coverage -abam #{bam_file.inspect} -b #{gff_file.inspect} -hist"
cmd1 += ' -s' unless options[:ignore_directions]
cov_lines_fwd = Bio::Commandeer.run cmd1, :log => log
#cov_lines_fwd = Bio::Commandeer.run "samtools view -b #{bam_file.inspect} 'gi|169887498|gb|CP000948.1|:1-2000' |bedtools coverage -abam - -b #{gff_file.inspect} -d -s", :log => log
if options[:ignore_directions]
  log.info "Parsing coverage profiles"
else
  log.info "Parsing forward aligned reads"
end
covs_fwd = get_covs.call(cov_lines_fwd)

unless options[:ignore_directions]
  cov_lines_rev = Bio::Commandeer.run "bedtools coverage -abam #{bam_file.inspect} -b #{gff_file.inspect} -hist -S", :log => log
  #cov_lines_rev = Bio::Commandeer.run "samtools view -b #{bam_file.inspect} 'gi|169887498|gb|CP000948.1|:1-2000' |bedtools coverage -abam - -b #{gff_file.inspect} -d -S", :log => log
  log.info "Parsing reverse aligned reads"
  covs_rev = get_covs.call(cov_lines_rev)
end

headers = [
  'contig',
  'type',
  'start',
  'end',
  'strand',
]
if options[:ignore_directions]
  headers.push 'average_coverage'
else
  headers.push 'forward_average_coverage'
  headers.push 'reverse_average_coverage'
end
headers.push 'annotation'
puts headers.join("\t")

covs_fwd.each do |feature, cov_fwd|
  cov_rev = covs_rev[feature] unless options[:ignore_directions]
  record = Bio::GFF::GFF3::Record.new(feature.join("\t"))

  products = record.attributes.select{|a| a[0] == 'product'}
  product = 'unannotated'
  if products.length == 1
    product = products[0][1]
  end

  to_print = [
    record.seqname,
    record.feature,
    record.start,
    record.end,
    record.strand,
    cov_fwd,
  ]
  to_print.push cov_rev unless options[:ignore_directions]
  to_print.push product
  puts to_print.join("\t")
end
