#!/usr/bin/env ruby

# Input a file output by
# fastacmd -T ...
# and output a file with lines only with
# sequence id\ttaxonomy id\tcommon name
#
#NCBI sequence id: gi|145225972|ref|YP_001136626.1|
#NCBI taxonomy id: 350054
#Common name: Mycobacterium gilvum PYR-GCK
#Scientific name: Mycobacterium gilvum PYR-GCK
#
seq_id = nil
taxon_id = nil
ARGF.each do |line|
  if matches = line.match(/NCBI sequence id\: (.*)/)
    seq_id = matches[1]
  end
  if matches = line.match(/NCBI taxonomy id\: (.*)/)
    taxon_id = matches[1]
  end
  if matches = line.match(/Common name\: (.*)/)
    name = matches[1]
    puts [
      seq_id, taxon_id, name
    ].join("\t")
  end
end