#!/usr/bin/env ruby


options = {
  :all_testing_data => '/home/ben/phd/plasmarithm/prediction_set_gathering/conservative_with_old_samples/training_set.50percent.manual.csv',
  :output_dir => '/home/ben/phd/plasmarithm/prediction_set_gathering/conservative_with_old_samples/cv'
}

Dir.mkdir options[:output_dir] unless File.exist?(options[:output_dir])
Dir.chdir options[:output_dir] do
  samples = File.open(options[:all_testing_data]).readlines
  
  outputs = (1..5).collect do |cv|
    File.open(File.join(options[:output_dir],"#{cv}.ids"),'w')
  end
  samples.each_with_index do |sample, index|
    output = index % 5
    outputs[output].print sample
  end
end
