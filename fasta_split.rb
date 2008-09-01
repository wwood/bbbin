#!/usr/bin/ruby 


# == Synopsis
#
# fasta_split: splits up a fasta file into several sections
#
# == Usage
#
# hello [OPTION] [FASTA_FILE]
#
# -h, --help:
#    show help
#
# --max-sequences x, -n x:
#    maximum number of sequences to be put in a single fasta file 
#    (NOT the number of files to be output)
#    default 10 000
#
# --prepend x, -p x:
#    add x to the start of the name of the file 
#    default 'split'
#
# --append x, -a x:
#    append x to the end of the file names
#    default 'Out'
#
# FASTA_FILE: The input fasta file (if not present, then STDIN is used).
#
# Example: fasta_split.rb --prepend 'number' --append 'file.fasta' -n2 <my.fasta
# This splits the my.fasta file into chunks of 2 sequences each, called 'number0000file.fasta', 'number0001file.fasta', etc.






require 'getoptlong'
require 'rdoc/usage'
require 'bio'

opts = GetoptLong.new(
		      [ '--help', '-h', GetoptLong::NO_ARGUMENT ],
		      [ '--max-sequences', '-n', GetoptLong::REQUIRED_ARGUMENT ],
		      [ '--prepend', '-p', GetoptLong::REQUIRED_ARGUMENT ],
		      [ '--append', '-a', GetoptLong::REQUIRED_ARGUMENT ]
		      )



# Parse arguments
prepend = 'split'
append = 'Out'
max_sequences = 10000
opts.each do |opt, arg|
  case opt
  when '--help'
    RDoc::usage
  when '--max-sequences'
    max_sequences = arg.to_i
  when '--prepend'
    prepend = arg
  when '--append'
    append = arg
  end
end


if ARGV.length > 1
  $stderr.puts "Failed to parse command line. More than 1 fasta file was specified"
  exit 1
end



# Read in the fasta file
if ARGV.length == 1
  #puts "Using file from the command line"
  @stream = File.open(ARGV.shift)
else
  #puts "Using stream"
  @stream = $stdin
end
fasta = Bio::FlatFile.open(Bio::FastaFormat, @stream)



# Print out the sequences
@entry = Bio::FastaFormat
file_number = 0
while @entry
  n = 0
  name = prepend+
    sprintf('%04d', file_number)+
    append
  file_number += 1


  output = Bio::FlatFile.new(Bio::FastaFormat, File.open(name, 'w'))

  while n<max_sequences and 
    n+=1
    @entry = fasta.next_entry

    break if !@entry

    output.to_io.puts @entry.to_s
  end
  
  output.close
end
