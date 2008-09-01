#!/usr/bin/ruby
require 'bio'



# Parse a FASTA file, trying to work out what is wrong with it

#TODO: when there is just a \r in a file, sequences are missed. I need to invent some way of detecting and fixing this. Will probably have to cache the sequence itself at the beginning, and do my own parsing before I give it to bioruby.



# Parse in file
stream = nil
if ARGV.length > 1
  $stderr.puts "only one sequence at a time sorry"
  exit 1
elsif ARGV.length == 1
  stream = File.open ARGV[0]
else
  stream = $stdin
end
fasta = Bio::FlatFile.new(Bio::FastaFormat, stream)


#hash of all the ids to check for duplicates
ids = {}
count = 0

s = fasta.next_entry

while s
   count += 1
  
  #if !(s.class==="Bio::FlatFile")
  #  puts "warning: strange 


  # Make sure the length is non-zero
  if s.length==0
    $stderr.puts "Zero length sequence: #{s.entry_id}"
  end


  if ids[s.entry_id]
    #TODO: this might actually generate duplicates, if there is one already created like that.
    newname = "#{s.entry_id}-#{ids[s.entry_id]}"
    ids[s.entry_id] += 1
    $stderr.puts "warning: duplicate entry id #{s.entry_id} changed to #{newname}"
    # s.entry_id = newname #why doesn't that work. Stupid.
    puts ">#{newname} #{s.comment}"
    puts "#{s.seq}"
  else
    ids[s.entry_id] = 1
    # Print out the parsed in sequence, in its new form
    puts s.to_s
  end

  
  

  s = fasta.next_entry
end

$stderr.puts "Counted #{count} sequences, including duplicates."
