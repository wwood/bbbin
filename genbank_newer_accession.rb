#!/usr/bin/env ruby

# Takes a file of genpept (full) / genbank format entries, and spits out what the newer gis, based on
# the 
# COMMENT     [WARNING] On Jan 12, 2007 this sequence was replaced by
#            gi:122065184.
#
# part of the paper. Spits it out in the form of accession tab replacing_entry
# and 

require 'rubygems'
require 'bio'

if __FILE__ == $0
  Bio::FlatFile.foreach(Bio::GenBank, ARGF) do |me|
    new_accession = nil
    
    comment_spaceless = me.comment.gsub(/\n\s*/,' ')
    if matches = comment_spaceless.match(/this sequence was replaced by (\S*)\./)
      new_accession = matches[1]
    end
    
    puts [
    me.accession,
    new_accession
    ].join("\t")
  end
end