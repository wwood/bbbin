#!/usr/bin/env ruby

require 'rubygems'
#require 'go' #provided by goruby
require 'eupathdb_gene_information_table' #provided by the reubypathdb gem
require '/home/ben/forays/goruby/lib/go'

if __FILE__ == $0
  inspected_go = 'GO:0005634'
  gor = Bio::Go.new
  
  # nucleus subsumer
  nucleus_subsumer = gor.subsume_tester(inspected_go)
  cc_subsumer = gor.subsume_tester('GO:0005575')
  # nucleus' ancestors or ancestors of any of its children. 
  # These annotations are ignored, because they aren't
  # positive or negative, merely 'cordial'
  exclusion_list = gor.cordial_cc(inspected_go)
  
  EuPathDBGeneInformationTable.new(ARGF).each do |g|
    go_terms = g.get_table('GO Terms').collect{|r| r['GO ID']}
    
    positive = false
    negative = false
    go_terms.each do |d|
      begin
        i = gor.primary_go_id d
        nuc = nucleus_subsumer.subsume?(i)
        cc = cc_subsumer.subsume?(i)
        
        positive = true if nuc 
        negative = true if cc and !(exclusion_list.include?(i)) and i!=inspected_go
      rescue RException
        # ignore when there are strange GO annotations, that don't seem to exist in the database
      end
    end
    
    # Print the result for this gene
    print "#{g['Gene Id']}\t"
    if positive
      if negative
        puts 'both'
      else
        puts 'positive'
      end
    elsif negative
      puts 'negative'
    else
      puts 'no data'
    end
  end
end