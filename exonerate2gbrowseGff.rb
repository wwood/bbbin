#!/usr/bin/env ruby
#
# Takes the output of an exonerate run (that has been run with the
# --showtargetgff or the other showgff and gives prints out just the GFF part of
# the file, formatted ready for gbrowse.
#
# Only print out the similarity part, and make sure each exon is printed
# differently

seqname = 'if this shows it is an error'

ARGF.each do |line|
  splits = line.split("\t");
  if splits.length == 9 and splits[2] == 'gene'
    stwos = splits[8].split(' ; ')
    seqname = stwos[1].match(/sequence (.+)/)[1]
  elsif splits.length == 9 and splits[2] == 'exon'
    puts [
      splits[0..1],
      'similarity',
      splits[3..7],
      "ID=#{seqname};Name=#{seqname}"
    ].flatten.join("\t")
  end
end