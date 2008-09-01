#!/usr/bin/ruby

require 'rubygems'
require 'bio'


#report = Bio::Blast::Report.new($stdin.readlines.join("\n"), :xml)

# for each hsp, print query name, hit name, starts and stops, e-value, and description of hit, separated by tabs
  puts [
    'Query Definition',
    'Hit Accession',
    'Hit Id',
    'Hit Definition',
    'HSP Query from',
    'HSP Query to',
    'HSP Hit from',
    'HSP Hit to',
    'HSP e-value'
  ].join("\t")

Bio::Blast.reports(ARGF) do |report|
  report.iterations.each do |iteration| 
    iteration.each {|hit|
      hit.each do |hsp|
        puts [
          hit.query_def ? hit.query_def : report.query_def,
          hit.accession,
          hit.hit_id,
          hit.definition,
          hsp.query_from,
          hsp.query_to,
          hsp.hit_from,
          hsp.hit_to,
          hsp.evalue
        ].join("\t")
      end
    }
  end
  #  report.each do |hit|
  #    print hit.target_id, "\t", hit.evalue, "\n" if hit.evalue < 0.001
  #  end
end


