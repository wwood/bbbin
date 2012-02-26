#!/usr/bin/env ruby

# GNU join confuses me, so I'm going to write something that doesn't require sorting, which is the usual case

separator = "\t"
join_field = 0

file1_hash = {}

# Gather file1 up
File.foreach(ARGV[0]) do |line|
  splits = line.chomp.split(separator)
  key = splits[join_field].strip
  if file1_hash[key] #assume keys are unique, advise and ignore subsequents otherwise
    $stderr.puts "Found duplicate key in file 1: `#{key}', ignoring after the first time"
    next
  else
    file1_hash[key] = splits
  end
end

# Go through file 2
file1_foundeds = {}
num_splits_in_file1 = nil
File.foreach(ARGV[1]) do |line|
  splits = line.chomp.split(separator)

  if num_splits_in_file1.nil?
    num_splits_in_file1 = splits.length #assume all lines have the same number of splits, fail otherwise
  elsif num_splits_in_file1 != splits.length
    raise Exception, "Unexpected number of splits over '#{separator}' found in the first file. I give up so I don't cause you pain later"
  end

  key = splits[join_field].strip
  file2_data = splits[1..splits.length-1]
  if file1_hash[key]
    if file1_foundeds[key].nil?
      puts [file1_hash[key],file2_data].flatten.join(separator)
      file1_foundeds[key] ||= 0
      file1_foundeds[key] += 1
    else
      $stderr.puts "Found key `#{key}' more than once in the second file. Fail."
    end
  else
    # Print empty cells for the first file
    (num_splits_in_file1).times {print separator}
    # Print the data for the second file
    puts file2_data.join(separator)
  end
end

# Print out all the remaining ones that aren't in file 2
file1_hash.each do |key, value|
  puts value.join(separator) unless file1_foundeds[key]
end
