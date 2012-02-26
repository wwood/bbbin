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
  num_splits_in_file1 ||= splits.length
  key = splits[join_field].strip
  if file1_hash[key]
    if file1_foundeds[key].nil?
      puts [file1_hash[key],splits[1..splits.length-1]].flatten.join(separator)
      file1_foundeds[key] ||= 0
      file1_foundeds[key] += 1
    else
      $stderr.puts "Found key `#{key}' more than once in the second file. Fail."
    end
  else
    $stderr.puts "Unable to find key `#{key}' in File 1 - I need file 1 to define all the keys in file 2, ignoring"
    next
  end
end

# Print out all the remaining ones that aren't in file 2
file1_hash.each do |key, value|
  puts value.join(separator) unless file1_foundeds[key]
end
