#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'

require 'bio-express_beta_diversity'
require 'descriptive_statistics'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
  :num_rarefactions => 100,
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <arguments>

Given an express beta diversity OTU table file, calculate one form of alpha diversity: the number of OTUs that are observed 2 or more times, in each sample\n\n"

  opts.on("-e", "--ebd FILE", "EBD file to read in [required]") do |arg|
    options[:ebd_file] = arg
  end
  opts.on("-c", "--count NUM",Integer, "Number of samples to do rarefaction with [required]") do |arg|
    options[:count] = arg
  end
  opts.separator "\nOptional parameters:\n\n"
  opts.on("--num-rarefactions NUM",Integer, "Number of times to do rarefaction [default: #{options[:num_rarefactions]}") do |arg|
    options[:num_rarefactions] = arg
  end

  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end; o.parse!
if ARGV.length != 0 or options[:ebd_file].nil? or options[:count].nil?
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)



otus = Bio::EBD::Format.parse_from_file(options[:ebd_file])
log.info "Read in #{otus.number_of_samples} sample for analysis"


# Get samples in a state for efficient-ish sampling
arrays_for_sampling = {}
otus.sample_counts.each do |sample_name, observations|
  sampled_array = []
  observations.each_with_index do |obs, i|
    otu = otus.otu_names[i]
    obs.to_i.times do
      sampled_array.push otu
    end
  end
  arrays_for_sampling[sample_name] = sampled_array
end


sample_names_to_observed_alphas = {}
options[:num_rarefactions].times do
  otus.sample_counts.each do |sample_name, its_obs|
    all_observations = arrays_for_sampling[sample_name].sample(options[:count])
    counts = {}
    all_observations.each do |otu_id|
      counts[otu_id] ||= 0
      counts[otu_id] += 1
    end
    alpha = counts.select{|otu_id, count| count >= 2}.length

    sample_names_to_observed_alphas[sample_name] ||= []
    sample_names_to_observed_alphas[sample_name].push alpha
  end
end

sample_names_to_observed_alphas.each do |sample_name, alphas|
  puts [
    sample_name,
    alphas.mean
    ].join("\t")
end