#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'csv'
require 'pp'

$:.unshift File.dirname(__FILE__)
require 'plasmarithm_common'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :input_data_directory => '/home/ben/phd/plasmarithm/prediction_set_gathering/conservative_with_old_samples/cv',
  :plasmo_int_groups_csv => File.join(File.dirname(__FILE__), 'input_data', 'nbt.1597-S3.csv'),
  :all_the_answers_are_at => '/home/ben/phd/plasmarithm/prediction_set_gathering/conservative_with_old_samples/training_set.50percent.manual.csv',
  :method => :whacky_unknown_thing, #This is a bad default, but eh.
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <arguments>

    Generate inputs that use the answers from other proteins. For the purposes of cross validation, some characterised proteins are off limits, see....\n"

  opts.on("-d", "--directory", "Directory where the cross validation data lives.
  Need 5 files (1.ids, 2.ids, 3.ids, 4.ids, 5.ids) in that directory that contain
  the PlasmoDB IDs in each of the 5 folds of cross-validation. [default: #{options[:input_data_directory]}]") do |f|
    options[:input_data_directory] = f
  end
  opts.on('-m', '--method METHOD', 'Use a different method for scoring each localisation for each PlasmoINT group') do |method|
    options[:method] = method.to_sym
  end
  opts.on('-o', '--output-dir DIR', 'Output the generated files to this directory') do |dir|
    raise unless File.exist?(dir) and File.directory?(dir)
    options[:output_data_directory] = dir
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

options[:output_data_directory] = options[:input_data_directory] if options[:output_data_directory].nil?

module Plasmarithm
  class CrossValidation
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
    def localisations
      %w(apical
      apicoplast
      cytosol
      endoplasmic_reticulum
      exported
      food_vacuole
      mitochondrion
      nucleus
      parasitophorous_vacuole
      plasma_membrane)
    end
    
    def log
      Bio::Log::LoggerPlus.new(LOG_NAME)
    end

    attr_accessor :partition_plasmodb_ids, :answers

    def read_plasmodb_ids_and_classifications(input_data_directory)
      # Read in the lists of IDs and their classification
      @partition_plasmodb_ids = {}
      @answers = {}
      i = 1
      while true
        filename = File.join(input_data_directory, "#{i}.ids")
        break unless File.exist?(filename)

        @partition_plasmodb_ids[i] = []
        CSV.foreach(filename, :col_sep => "\t") do |row|
          plasmodb = row[0]
          answer = row[1].gsub(' ','_')
          @partition_plasmodb_ids[i].push plasmodb
          raise unless answers[plasmodb].nil?
          @answers[plasmodb] = answer
        end
        log.debug "Found #{@partition_plasmodb_ids[i].length} PlasmoDB IDs in cross validation IDs file #{filename}"

        i += 1
      end
    end

    def check_cross_validation_groups
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
    end

    attr_accessor :plasmoint_groups_to_ids, :plasmoint_ids_to_groups

    def read_plasmoint_data(plasmoint_groups_file)
      # Read in the PlasmoInt data. IDs => group number hash, group number => IDs hash
      @plasmoint_groups_to_ids = {}
      @plasmoint_ids_to_groups = {}
      line_no = 1
      CSV.foreach(plasmoint_groups_file, :col_sep => "\t") do |row|
      # First 2 lines are headers
        unless line_no > 2
        line_no += 1
        next
        end

        group_id = row[0].to_i
        raise unless group_id > 0
        plasmodb = row[1]

        @plasmoint_groups_to_ids[group_id] ||= []
        @plasmoint_groups_to_ids[group_id].push plasmodb

        @plasmoint_ids_to_groups[plasmodb] = group_id
      end
    end

    # A representation of a training and testing set tuple
    class CrossSection
      def initialize(cross_validation_object, training_samples, testing_samples)
        @cv = cross_validation_object
        @training_samples = training_samples
        @testing_samples = testing_samples
        @ranking_method = :whacky_unknown_thing
      end
      
      def log
        Bio::Log::LoggerPlus.new(LOG_NAME)
      end
      
      attr_accessor :ranking_method

      attr_reader :loc_background_counts, :loc_background_percents

      def compute_background_percents
        @loc_background_counts = {}
        @training_samples.each do |plasmodb|
          if @cv.answers[plasmodb].nil?
            raise "I think #{plasmodb} should have known localisation, but it doesn't seem to.. fail"
          end
          @loc_background_counts[@cv.answers[plasmodb]] ||= 0
          @loc_background_counts[@cv.answers[plasmodb]] += 1
        end
        @loc_background_percents = {}
        @loc_background_counts.each do |loc, count|
          @loc_background_percents[loc] = count.to_f/@loc_background_counts.values.sum
        end
        log.info "Found background percentages #{loc_background_percents}"
      end

      # Print a particular sample to the output IO given
      def print_sample(plasmodb, output)
        group_id = @cv.plasmoint_ids_to_groups[plasmodb]

        group_loc_counts, unknown_count, plasmoint_id = gather_plasmoint_data(plasmodb)
        log.debug "Known classifications for this group: #{group_loc_counts}"

        # print headers
        output.print plasmodb

        # for each of the localisations, output the percentage
        @cv.localisations.each do |loc|
        # output the percent in this group that have that localisation
          output.print "\t"
          if @ranking_method == :whacky_unknown_thing
            output.print whacky_unknown_thing(loc, group_loc_counts, unknown_count, plasmoint_id)
          elsif @ranking_method == :simple_maximal
            output.print simple_maximal_rank(loc, group_loc_counts, unknown_count, plasmoint_id)
          else
            raise "unknown ranking method: #{@ranking_method}"
          end
        end
        output.print "\t"
        output.print @cv.answers[plasmodb]
        output.puts
      end

      # Given a plasmodb identifier, return percentages of the group that contains the given plasmodb identifier
      def gather_plasmoint_data(plasmodb)
        group_loc_counts = {}
        @cv.localisations.each{|loc| group_loc_counts[loc]=0}
        plasmoint_id = @cv.plasmoint_ids_to_groups[plasmodb]

        unknown_count = nil
        unless plasmoint_id.nil?
          @cv.plasmoint_groups_to_ids[plasmoint_id].each do |group_plasmodb|
            next if group_plasmodb == plasmodb # Don't count the training sample itself
            answer = @cv.answers[group_plasmodb]
            if answer
            group_loc_counts[answer] += 1
            end
          end
          log.debug "Numbers of members of the group: #{@cv.plasmoint_groups_to_ids[plasmoint_id].length}"
          # unknowns = total in group - knowns - this_sample
          unknown_count = @cv.plasmoint_groups_to_ids[plasmoint_id].length - group_loc_counts.values.sum - 1
          log.debug "Found unknown count #{unknown_count}"
        end

        return group_loc_counts, unknown_count, plasmoint_id
      end
      
      # Simply return the number of other proteins in the group that have the given localisation.
      # To predict based on this, one (perhaps) simply chooses that localisation that is highest
      def simple_maximal_rank(loc, group_loc_counts, unknown_count, plasmoint_id)
        if plasmoint_id.nil? or group_loc_counts.values.sum == 0
          # If there is no group, or the group has no characterised members 
          return loc_background_percents[loc]
        else
          # Some thing known. Report what is known
          return group_loc_counts[loc]
        end
      end

      # Return the final output for a given localisation and plasmodb. Perhaps overly context-dependent, but oh well.
      def whacky_unknown_thing(loc, group_loc_counts, unknown_count, plasmoint_id)
        if plasmoint_id.nil?
          return loc_background_percents[loc]
        else
          total = @cv.plasmoint_groups_to_ids[plasmoint_id].length
          num_this_localisation = loc_background_percents[loc]*unknown_count + group_loc_counts[loc]
          return num_this_localisation.to_f/total
        end
      end
    end
  end
end

cv = Plasmarithm::CrossValidation.new
cv.read_plasmodb_ids_and_classifications options[:input_data_directory]
cv.check_cross_validation_groups
cv.read_plasmoint_data options[:plasmo_int_groups_csv]

# For each of the ID lists
cv.partition_plasmodb_ids.each do |testing_partition_id, testing_plasmodb_ids|
# This iterated set is the testing set. Use the rest of the sets to train.
# output to trainingX.csv
  output = File.open(File.join(options[:output_data_directory], "testing#{testing_partition_id}.training_data.csv"), 'w')
  #output = $stderr

  # Output headers
  headers = [
    'plasmodb',
    cv.localisations,
    'answer',
  ].flatten
  output.puts headers.join("\t")

  # work out the background prevalence of each localisation in this training set
  training_sets = cv.partition_plasmodb_ids.keys.reject{|t| t==testing_partition_id}
  log.debug "Creating traning set, where the testing set is #{testing_partition_id} - training sets #{training_sets.inspect}"

  training_samples = training_sets.collect{|t_id| cv.partition_plasmodb_ids[t_id]}.flatten
  log.debug "Using #{training_samples.length} training samples for the crossval that tests set #{testing_partition_id}"

  testing_samples = cv.partition_plasmodb_ids[testing_partition_id]
  cross_section = Plasmarithm::CrossValidation::CrossSection.new(cv, training_samples, testing_samples)
  cross_section.compute_background_percents
  
  cross_section.ranking_method = options[:method]

  # for each of the ids in this set,
  training_samples.each do |plasmodb|
    cross_section.print_sample(plasmodb, output)
  end

  output.close

  # Print the testing sample data, which is only based on the data extracted from
  File.open(File.join(options[:output_data_directory], "testing#{testing_partition_id}.testing_data.csv"), 'w') do |testing_output|
    testing_output.puts headers.join("\t")

    testing_samples = cv.partition_plasmodb_ids[testing_partition_id]
    testing_samples.each do |plasmodb|
      cross_section.print_sample(plasmodb, testing_output)
    end
    
    testing_output.close
  end
end