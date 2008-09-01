#!/usr/bin/ruby

# Convert a csv file into a FASTA file, where a single column of the FASTA is the id,
# and another column is the sequence

require 'csv'
require 'optparse'
require 'ostruct'

id_col = 5
sequence_col = 14




options = OpenStruct.new
OptionParser.new do |opts|
  opts.banner = "Usage: csv2fasta.rb [options]"
  
  opts.on("-i", "--id-column [INTEGER]", "Specify Column in CSV corresponding to clone id. Default 5, first is 0.") do |v|
    id_col = v.to_i
  end

  opts.on("-s", "--sequence-column [INTEGER]","Specify Column in CSV corresponding to sequence. Default 14, first is 0.") do |v|
    sequence_col = v.to_i
  end
end.parse!

first = true
count = 0
hash = Hash.new

CSV::Reader.parse(File.open(ARGV[0], 'rb')) do |row|
  
  count += 1
  
  # skip the first row, and blank lines
  if count == 1 || !row ||row==='' ||row.length==0
    next
  end
  
  id = row[id_col]
  sequence = row[sequence_col]
  
  # Record the id and sequence, and compare to the last one if any exists
  blank = 'BLANK!!!'
  if !hash[id] #no entry yet
    # create entry
    if !sequence || sequence === ''
      hash[id] = blank
    else
      hash[id] = sequence
    end
    
  elsif hash[id]===blank #blank entry already
    if sequence
      $stderr.puts "WARNING: line #{count}: Blank and non-blank sequences define the same clone: #{id}"
    end
    
  else #sequence already
    # Ensure the sequences are the same
    if hash[id] === sequence
      next
    else
      $stderr.puts "WARNING: line #{count}: Duplicate names with different sequences: #{id}. Sequences '#{hash[id]}' and '#{sequence}'"
      next
    end
  end
  
  #$stderr.puts "ROW: #{row}"
  #$stderr.puts
  
  if !id
    $stderr.puts "WARNING: line #{count}: Nil ID"
  end
  
  if sequence
    match = id.match(/^Ascid Clone (\d+) *$/)
    if match
      puts ">#{id.gsub(/ /,'_')}"
      puts sequence
    else
      $stderr.puts "WARNING: line #{count}: Malformed ID: '#{id}'"
    end
  end
  
  
end
