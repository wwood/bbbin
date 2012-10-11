#!/usr/bin/env ruby

# Take a fasta file and list of blast hits, and annotate the FASTA by using the best blast hit method.

require 'bio'
require 'optparse'
require 'open3'

options = {
  :blast_db_path => '/srv/whitlam/bio/db/ncbi/nr', #where to get the sequence information from
  :min_evalue => 1e-5,
}
OptionParser.new do |opts|
  opts.banner = "Usage: annotate_fasta.rb -f <fasta_path> -b <blast_results_path>"

  opts.on('-f', "--fasta FASTA_TO_ANNOTATE", "Fasta file that needs annotating") do |v|
    options[:fasta] = v
  end

  opts.on('-b', "--blast BLAST_FILE", "BLAST of fasta against the db (-outfmt 6)") do |v|
    options[:blast] = v
  end

  opts.on('-s', "--species SPECIES_NAME", "Name of the species/community/genome build e.g. 'RCIIv2.3.1_'. You probably want an underscore at the end so that the gene IDs are separated from the species name visually, yet still both are before the first space.") do |v|
    options[:species_name] = v
  end

  opts.on('-d', "--blast-db PATH", "Path to the blast DB where the information is gathered [default: #{options[:blast_db_path]}]") do |v|
    options[:blast_db_path] = v
  end
end.parse!

raise Exception, "Please specify a species/community name (-s/--species)" unless options[:species_name]
raise unless options[:fasta]
raise unless options[:blast]

query_to_hit = {}

# cache the blast hits
File.open(options[:blast]).each_line do |line|
  splits = line.chomp.split("\t")

  query = splits[0]
  next if query_to_hit[query]

  hit = splits[1]
  evalue = splits[10]
  next if evalue.to_f > options[:min_evalue]

  query_to_hit[query] = hit
end

# Run through the fasta file, annotating as you go along
Bio::FlatFile.foreach(options[:fasta]) do |s|
  name = 'hypothetical'
  if query_to_hit[s.entry_id]
    # extract best blast hit full name
    cmd = "blastdbcmd -entry '#{query_to_hit[s.entry_id]}' -db #{options[:blast_db_path]} -outfmt '%t'"

    # execute command and capture both stdout, and stderr
    Open3.popen3(cmd) do |stdin, stdout, stderr|
      error = stderr.readlines
      if error.length > 0
        raise Exception, "Error while extracting info from blastdb:\n#{error}"
      end

      name = stdout.readlines[0].gsub(/\[.*/,'')

    end
  end
  print ">#{options[:species_name]}#{s.entry_id} "
  puts name
  puts s.seq
end

