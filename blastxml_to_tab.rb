#!/usr/bin/env ruby

require 'rubygems'
require 'bio'
require 'optparse'

USAGE = "Usage: blastxml_to_tab.rb [-n] <blastxml_file>"
options = {
  :ncbi => false
}
OptionParser.new do |opts|
  opts.banner = USAGE

  opts.on("-n", "--ncbi", "Print out tabular like NCBI (incomplete at the moment)") do |v|
    options[:ncbi] = v
  end
end.parse!

#report = Bio::Blast::Report.new($stdin.readlines.join("\n"), :xml)

# for each hsp, print query name, hit name, starts and stops, e-value, and description of hit, separated by tabs
if options[:ncbi]
  puts [
    '# Query ID',
    'Subject ID',
    '% identity',
    'alignment length',
    'mismatches',
    'gap openings',
    'HSP Query from',
    'HSP Query to',
    'HSP Hit from',
    'HSP Hit to',
    'HSP e-value'
  ].join("\t")
else
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
end
Bio::Blast.reports(ARGF) do |report|
  report.iterations.each do |iteration| 
    iteration.each {|hit|
      hit.each do |hsp|
        if options[:ncbi]
          puts [
            hit.query_def ? hit.query_def : report.query_def,
            hit.hit_id,
            hsp.identity,
            hsp.align_len,
            hsp.mismatch_count,
            hsp.gaps,
            hsp.query_from,
            hsp.query_to,
            hsp.hit_from,
            hsp.hit_to,
            hsp.evalue
          ].join("\t")
        else
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
      end
    }
  end
  #  report.each do |hit|
  #    print hit.target_id, "\t", hit.evalue, "\n" if hit.evalue < 0.001
  #  end
end


