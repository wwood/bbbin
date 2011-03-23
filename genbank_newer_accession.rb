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

# monkey in the comment to new gi method
class Bio::GenBank
  # Return the newer GI for this entry, or the current gi if there are no newer entries. This is encoded in the comment section, using "this sequence was replaced by gi:\d+" where the \d+ is returned.
  # options:
  # recursive: contact NCBI using the efetch mechanism, and continue getting the more up to date accession until no more updates are found.
  # verbose: update when new GIs are found 
  def updated_gi(options = {})
    options ||= []
    
    update_gi_from_local_entry = lambda do |gb|
      # remove the breaking of the comment over several lines
      comment_spaceless = gb.comment.gsub(/\n\s*/,' ')
      
      if matches = comment_spaceless.match(/this sequence was replaced by (\S*)\./)
        new_accession = matches[1]
        $stderr.puts "Found an updated GI: #{new_accession}"
        
        # go through it recursively. Is that newer identifier still out of date?
        if matches2 = matches[1].gsub(/gi:/,'')
          if options[:verbose]
            $stderr.puts "Found an updated GI: #{matches2}"
          end
          return matches2
        end
        
        # parsing error of comment
        raise ParseException, "Could not parse comment `#{matches[1]}' from the comment section of #{accession}: #{comment}"
      else
        # no updated entry
        return gb.gi
      end
    end # update_gi_from_local_entry lambda
    
    if options[:verbose]
      $stderr.puts "Determining updated GIs for GenBank accession #{self.accession}"
    end
    
    newer_gi = update_gi_from_local_entry.call(self)
    if options[:recursive]
      puts self.class
      gb = self
      # while there are newer gi's being found
      while newer_gi != gb.gi
        gb = Bio::GenBank.new(Bio::NCBI::REST.efetch(newer_gi, {:rettype => 'gb', :db => ENTREZ_DATABASE})) # is there no way to figure this out automatically?
        newer_gi =  update_gi_from_local_entry.call gb
      end
    end
    return newer_gi
  end
end




if __FILE__ == $0
  Bio::NCBI.default_email = "go@away.com"
  
  Bio::FlatFile.foreach(Bio::GenBank, ARGF) do |me|
    new_accession = nil
    
    updated_gi = me.updated_gi({:recursive => true, :verbose => :true})
    if updated_gi != me.gi
      new_accession = updated_gi
    end
    
    puts [
    me.accession,
    new_accession
    ].join("\t")
  end
end