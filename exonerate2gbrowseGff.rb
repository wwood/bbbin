#!/usr/bin/env ruby
#
# Takes the output of an exonerate run (that has been run with the
# --showtargetgff or the other showgff and gives prints out just the GFF part of
# the file, formatted ready for gbrowse.
#
# Only print out the similarity part, and make sure each exon is printed
# differently

ARGF.each do |line|
  splits = line.split("\t");
  if splits.length == 9 and splits[2] == 'similarity'
    stwos = splits[8].split(' ; ')
    seqname = stwos[1].match(/Query (.+)/)[1]
    stwos[2..(stwos.length-1)].each do |align|
      matches = align.match(/Align (\d+) \d+ (\d+)/)
      raise unless matches
      pilon = matches[1].to_i
      length = matches[2].to_i
      start = 0
      stop = 0
      if splits[6] == '-'
        start = pilon-length
        stop = pilon+1
      elsif splits[6] == '+'
        start = pilon
        stop = pilon+length-1
      end

      puts [
        splits[0..2],
        start,
        stop,
        splits[5..7],
        "ID=#{seqname};Name=#{seqname}"
      ].flatten.join("\t")
    end
  end
end