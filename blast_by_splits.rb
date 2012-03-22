#!/usr/bin/env ruby

# Given a fasta file, blast it against a database. However, first split the query sequences up equally so that the load can be spread over multiple CPUs

require 'rubygems'
require 'tempfile'
require 'bio'
require 'optparse'


options = {
  :threads => 24,
}
OptionParser.new do |opts|
  opts.banner = "Usage: blast_by_splits.rb --query <query_fasta> --db <blast_database_path>"

  opts.on('-i', "--query QUERY", "Query fasta file") do |v|
    options[:query] = v
  end
  
  opts.on('-d', "--db DB", "path to blast database") do |v|
    options[:db] = v
  end

  opts.on('-a', "--threads NUM_THREADS", "Number of CPUs to spread the load across") do |v|
    options[:threads] = v.to_i
  end
end.parse!

raise unless options[:query]
raise unless options[:db]
raise unless options[:threads] > 0

input_temps = (1..options[:threads]).collect{Tempfile.new('blast_by_splitsIn')}


# Spread input sequences across thread files
total_sequences = 0
Bio::FlatFile.open(options[:query]).each_with_index do |s, i|
  total_sequences += 1
  input_temps[i % options[:threads]].puts s.to_s
end
input_temps.each do |temp|
  temp.close
end
$stderr.puts "Starting BLASTs of #{total_sequences} sequences across #{options[:threads]} CPUs .."
output_temps = (1..options[:threads]).collect{Tempfile.new('blast_by_splitsOut')}

# Start each of the threads
blast_threads = (1..options[:threads]).collect do |i|
  Thread.new do
    `blastp -query '#{input_temps[i-1].path}' -db '#{options[:db]}' -num_descriptions 1 -num_alignments 1 -outfmt 6 -out #{output_temps[i-1].path}`
  end
end

#wait for the all to be finished, then print all the results out as they finish
blast_threads.each_with_index do |thread, i|
  thread.join
    File.open(output_temps[i]).each_line do |line| 
    print line
  end
end

