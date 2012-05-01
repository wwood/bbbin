#!/usr/bin/env ruby

# Given a fasta file, blast it against a database. However, first split the query sequences up equally so that the load can be spread over multiple CPUs

require 'rubygems'
require 'tempfile'
require 'bio'
require 'optparse'


options = {
  :threads => 24,
  :max_target_seqs => 1,
  :outfmt => 6,
}
o = OptionParser.new do |opts|
  opts.banner = "Usage: blast_by_splits.rb --query <query_fasta> --db <blast_database_path> --out <output_filename>"

  opts.on('-i', "--query QUERY", "Query fasta file [required]") do |v|
    options[:query] = v
  end
  
  opts.on('-d', "--db DB", "path to blast database [required]") do |v|
    options[:db] = v
  end
  
  opts.on('-o', "--out OUTPUT_FILE", "File to dump blast results to. Note that the sequences with blast hits are probably not in the same order as the input file [required]") do |v|
    options[:output_file] = v
  end

  opts.on('-a', "--threads NUM_THREADS", "Number of CPUs to spread the load across [default: #{options[:threads]}]. Note that this doesn't use BLAST's in-built option, but splits up the query file into chunks and uses 1 cpu per chunk (and 1 call to BLAST per chunk, and therefore the db is loaded into memory by BLAST once for each chunk)") do |v|
    raise Exception, "Number of threads (-a/--threads) must be > 0, quitting" unless options[:threads] > 0
    options[:threads] = v.to_i
  end
    
  opts.on('-t', "--max_target_seqs NUM_SEQS", "How many hits to report for each query sequence? [default: #{options[:max_target_seqs]}]") do |v|
    options[:max_target_seqs] = v.to_i
  end

  opts.on('-m', "--outfmt FORMAT_NUMBER", "Output format [default: #{options[:max_target_seqs]} (tab-separated values)]") do |v|
    options[:outfmt] = v
  end

  opts.on('-e', "--evalue E_VALUE", "E-value cutoff [default: use the default in BLAST (10, I think)]}") do |v|
    options[:evalue] = v
  end
end.parse!

unless options[:query] and options[:db] and options[:output_file]
  $stderr.puts "I need -i/--query, -d/--db and -o/--out arguments before I can run, quitting\n\n"
  $stderr.puts o.usage
  exit 1
end

raise Exception, "Unexpected unqualified arguments: #{ARGV.join(' ')}" unless ARGV.length == 0

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
    num_sequences_argument = "-max_target_seqs #{options[:max_target_seqs]}"
    if %w(0 1 2 3 4 5).include?(options[:outfmt].to_s)
      num_sequences_argument = "-num_descriptions #{options[:max_target_seqs]} -num_alignments #{options[:max_target_seqs]}"
    end
    cmd = "blastp -query '#{input_temps[i-1].path}' -db '#{options[:db]}' #{num_sequences_argument} -outfmt #{options[:outfmt]} -out #{output_temps[i-1].path}"
    cmd = "#{cmd} -evalue #{options[:evalue]}" unless options[:evalue].nil?
    #$stderr.puts "Running: #{cmd}"
    `#{cmd}`
  end
end

# Setup the output file, and try writing to it so it fails asap, so you don't have to wait for the blast to finish to find out there's a problem
out = File.open(options[:output_file],'w')
out.print ''

#wait for the all to be finished, then print all the results out as they finish
blast_threads.each_with_index do |thread, i|  
  thread.join
    File.open(output_temps[i]).each_line do |line| 
    out.print line
  end
end

out.close

