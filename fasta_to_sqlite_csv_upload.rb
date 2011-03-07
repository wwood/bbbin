#!/usr/bin/env ruby

# Takes a fasta file as input, and then output CSV files so that it can be uploaded using the SQlite (and perhaps other DBs) csv import mechanism. Output files:
# 1. CSV file for names of each sequence (id, sequence_source_id, name)
# 2. CSV file for amino acid sequences themselves (id, sequence_name_id, sequence)
# 3. code for upload of the csv to be input into SQLite's command line (sqlite specific)

require 'bio'

unless ARGV.length > 0 and ARGV.length%2 == 0
  $stderr.puts "Usage: fasta_to_sqlite_csv_upload.rb [<sequence_source_description> <fasta_file> ..] <sequence_source_description> <fasta_file>"
  exit 1
end

# global variables that may be subject to option parsing
starting_sequence_name_id = 1
starting_sequence_id = 1 #for the sequence itself, not the name of the sequence
sequence_source_id = 1
upload_date = "2011-02-28 13:32:37"
sequence_names_csv_filename = 'sequence_names_to_upload.csv'
sequence_csv_filename = 'sequence_to_upload.csv'
separator = "\t"

# open CSV files to write to
name_csv_file = File.open(sequence_names_csv_filename,'w')
sequence_csv_file = File.open(sequence_csv_filename,'w')



# A re-usable method for the upload of a fasta file-source_name combination, so that multiple fasta files can be uploaded in order
upload_fasta = lambda do |source_name,  fasta_filename|
  # insert the sequence source
  puts "insert into sequence_sources values(#{sequence_source_id}, '#{source_name}','#{upload_date}','#{upload_date}');"
                                          
                                          uploaded = 0
                                          Bio::FlatFile.foreach(fasta_filename) do |seq|
    # print name of sequence
    name_csv_file.puts [
    starting_sequence_name_id,
    seq.entry_id,
    sequence_source_id,
    upload_date, upload_date
    ].join(separator)
    
    
    # print sequence itself
    sequence_csv_file.puts [
    starting_sequence_id,
    starting_sequence_name_id,
    seq.seq,
    upload_date, upload_date
    ].join(separator)
    
    # increment for end of loop, for next time
    starting_sequence_name_id += 1
    starting_sequence_id += 1
    
    uploaded += 1
  end
  sequence_source_id += 1 #for the next source
  $stderr.puts "Converted #{uploaded} sequences to CSV format for upload from #{source_name}, #{fasta_filename}"
end



# One pair of fasta files at once
while ARGV.length > 0
  current_fasta_filename = ARGV.pop
  current_source_name = ARGV.pop
  
  upload_fasta.call(current_source_name, current_fasta_filename)  
end

# Only need to import once for the whole show

# define the separator for the CSV upload
puts '.separator "\t"'
# upload the sequence names
puts ".import #{sequence_names_csv_filename} coding_regions"
# upload the sequences
puts ".import #{sequence_csv_filename} amino_acid_sequences"
