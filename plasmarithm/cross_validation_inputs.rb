#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'csv'
require 'plasmarithm_common'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :input_data_directory => nil,
  :plasmo_int_groups_csv => File.join('input_data', 'nbt.1597-S3.csv'),
  :all_the_answers_are_at => '/home/ben/phd/plasmarithm/prediction_set_gathering/conservative_with_old_samples/training_set.50percent.manual.csv'
}
o = OptionParser.new do |opts|
#TODO Fill in usage, description and option parsing below
  opts.banner = "
    Usage: #{SCRIPT_NAME} <arguments>

    Generate inputs that use the answers from other proteins. For the purposes of cross validation, some characterised proteins are off limits, see....\n"

  opts.on("-d", "--directory", "Directory where the cross validation data lives.
  Need 5 files (1.ids, 2.ids, 3.ids, 4.ids, 5.ids) in that directory that contain
  the PlasmoDB IDs in each of the 5 folds of cross-validation. [REQUIRED]") do |f|
    options[:operation] = OVERALL
  end

  # logger options
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") do |q|
    Bio::Log::CLI.trace('error')
  end
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") do | name |
    options[:logger] = name
  end
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG") do | s |
    Bio::Log::CLI.trace(s)
  end
end
o.parse!
if ARGV.length != 0
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]) #bio-logger defaults to STDERR not STDOUT, I disagree
log = Bio::Log::LoggerPlus.new(LOG_NAME)
Bio::Log::CLI.configure(LOG_NAME)

options[:output_data_directory] = options[:input_data_directory]

# Setup column names
# ben@uyen:20120602:~/bin/plasmarithm$ awk -F'    ' '{print $2}' ~/phd/plasmarithm/prediction_set_gathering/conservative_with_old_samples/training_set.50percent.manual.csv |sort |uniq -c
# 44 apical
# 22 apicoplast
# 23 cytosol
# 17 endoplasmic_reticulum
# 30 exported
# 3 food vacuole
# 15 mitochondrion
# 22 nucleus
# 4 parasitophorous vacuole
# 6 plasma membrane

localisations = %w(apical
apicoplast
cytosol
endoplasmic_reticulum
exported
food_vacuole
mitochondrion
nucleus
parasitophorous_vacuole
plasma_membrane)

answers = {}
CSV.foreach(options[:all_the_answers_are_at],col_sep => "\t") do |row|
  answer = row[1].gsub(' ','_')
  raise unless localisations.include? answer
  answers[row[0]] = answer
end

# Read in the lists of IDs
partition_plasmodb_ids = {}
i = 1
while true
  filename = File.join(options[:input_data_directory], "#{i}.ids")
  break unless File.exist?(filename)

  partition_plasmodb_ids[i] = File.open(filename).read.split(/\s+/)
  log.debug "Found #{partition_plasmodb_ids[i]} PlasmoDB IDs in cross validation IDs file #{filename}"

  i += 1
end

lengths = partition_plasmodb_ids.collect{|num, ids| ids.length}
min = lengths.min
max = lengths.max
log.info "Found #{partition_plasmodb_ids.length} cross validation groups, with min #{min} and max #{max} identifiers in them"

# check that they are all unique ids
unless partition_plasmodb_ids.values.flatten.length == partition_plasmodb_ids.values.flatten.uniq.length
  raise "Found some identifiers that were specified multiple times"
end
# check that they all have pretty much the same number of ids in them
unless min==max or min==max-1
  raise "Unexpected numbers of identifiers found in some group"
end

# Read in the PlasmoInt data. IDs => group number hash, group number => IDs hash
plasmoint_groups_to_ids = {}
plasmoint_ids_to_groups = {}
line_no = 1
CSV.foreach(File.open(options[:plasmo_int_groups_csv]), :col_sep => "\t") do |row|
  # First 2 lines are headers
  unless line_no > 2
    line_no += 1
    next
  end
  
  group_id = row[0].to_i
  raise unless group_id > 0
  plasmodb = row[1]

  plasmoint_groups_to_ids[group_id] ||= []
  plasmoint_groups_to_ids[group_id].push plasmodb

  plasmoint_ids_to_groups[plasmodb] = group_id
end

# For each of the ID lists
partition_plasmodb_ids.each do |testing_partition_id, testing_plasmodb_ids|
# This iterated set is the testing set. Use the rest of the sets to train.
# output to trainingX.csv
  output = File.open(File.join(options[:output_data_directory], "testing#{testing_partition_id}.training_data.csv"))

  # Output headers
  output.puts [
    'plasmodb',
    localisations
  ].flatten.join(',')

  # work out the background prevalence of each localisation in this training set
  training_sets = partition_plasmodb_ids.keys.reject{|t| t==testing_partiting_id}
  training_samples = partition_plasmodb_ids[training_sets].flatten
  log.debug "Using #{training_samples.length} training samples for the crossval that tests set #{testing_partition_id}"
  loc_background_counts = {}
  training_samples.each do |plasmodb|
    if answers[plasmodb].nil?
      raise "I think #{plasmodb} should have known localisation, but it doesn't seem to.. fail"
    end
    loc_background_counts[answers[plasmodb]] ||= 0
    loc_background_counts[answers[plasmodb]] += 1
  end
  loc_background_percents = loc_background_counts.collect{|c| c.to_f/loc_background_counts.sum}

  # for each of the ids in this set,
  training_samples.each do |plasmodb|
    print plasmodb
      
    group_id = plasmoint_ids_to_groups[plasmodb]

    group_loc_counts = {}
    localisations.each{|loc| group_loc_counts[loc]=0}
    plasmoint_groups_to_ids[plasmodb].each do |group_plasmodb|
      next if group_plasmodb == plasmodb # Don't count the training sample itself
      answer = answers[group_plasmodb] 
      if answer
        group_loc_counts[answer] += 1
      end
    end
    # unknowns = total in group - knowns - this_sample
    unknown_count = plasmoint_groups_to_ids[plasmodb].length - group_loc_counts.values.sum - 1
    
    # for each of the localisations, output the percentage
    localisations.each do |loc|
      # output the percent in this group that have that localisation
      total = plasmoint_groups_to_ids[plasmodb].length
      
      print "\t"
      print loc_background_percents[loc]*unknown_count + group_loc_counts[loc]
    end
    puts
  end

  output.close
end