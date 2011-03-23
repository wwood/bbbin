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

ENTREZ_DATABASE = 'protein'
Bio::NCBI.default_email = "go@away.com"

# monkey in the comment to new gi method
class Bio::GenBank
  # Return the newer GI for this entry, or nil if none exist. This is encoded in the comment section, using "this sequence was replaced by gi:\d+" where the \d+ is returned. 
  def updated_gi
    # remove the breaking of the comment over several lines
    comment_spaceless = me.comment.gsub(/\n\s*/,' ')
    
    if matches = comment_spaceless.match(/this sequence was replaced by (\S*)\./)
      new_accession = matches[1]
      
      # go through it recursively. Is that newer identifier still out of date?
      if matches2 = matches[1].gsub(/gi:/,'')
        return newer_gi
      end
      
      # parsing error of comment
      raise ParseException, "Could not parse comment `#{matches[1]}' from the comment section of #{accession}: #{comment}"
    end
  end
end




if __FILE__ == $0
  Bio::FlatFile.foreach(Bio::GenBank, ARGF) do |me|
    new_accession = me.updated_gi
    
    if new_accession
      newer_returned = Bio::NCBI::REST.efetch(newer_gi, {:rettype => 'gb', :db => ENTREZ_DATABASE})
    end
    
    puts [
    me.accession,
    new_accession
    ].join("\t")
  end
end