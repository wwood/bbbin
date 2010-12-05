#!/usr/bin/env ruby

# Takes in a tab separated blast output file, and returns a 3 column
# tab separated file for use with sequenceChop.pl
#
# The input blast file may contains hits from a number of query proteins
# to a number of hits. This file collates the hits in 2 ways:
# * One-way blastclust - hit and query sequences are lumped together if they hit the same thing
# * hit extension - if 2 query sequences hit the same domain but with different lengths, the earliest start and latest stop will be included 

require 'rubygems'
require '/home/ben/forays/domainFinding/1/domain_finder_objects'
require 'fastercsv'
require 'optparse'

options = ARGV.getopts('q')

USAGE = [
  'Usage: blast_to_chopping.rb [-q] <blast_result_csv_file>',
  "\t-q: print domains from the (presumed blastx) query, merging overlapping domains except when their frames are different. Default is to print the hit sequences only."
].join("\n")

options['q'] ||= false #by default, print hits, not queries

unless ARGV.length == 1
  $stderr.puts USAGE
  exit
end

if __FILE__ == $0
  # An array of arrays
  clusters = []
  domains = {} #hit name to protein from domain_finder_objects
  
  # Read in blast file
  FasterCSV.foreach(ARGV[0], :col_sep => "\t") do |row|
    query = row[0]
    hit = row[1]
    query_start = row[6].to_i
    query_stop = row[7].to_i
    hit_start = row[8].to_i
    hit_stop = row[9].to_i
    
    # For each line, do the one-way blast clustering by merging the current hit into the cluster
    # search through the clusters looking for hitting query and hit
    query_cluster_index = nil
    clusters.each_with_index do |cluster, i|
      if cluster.include?(query)
        raise Exception, "Programming error - found 2 clusters with #{query}" unless query_cluster_index.nil?
        query_cluster_index = i
      end
    end
    hit_cluster_index = nil
    clusters.each_with_index do |cluster, i|
      if cluster.include?(hit)
        raise Exception, "Programming error - found 2 clusters with #{hit}" unless hit_cluster_index.nil?
        hit_cluster_index = i
      end
    end
    
    # Merge / Create clusters
    if !(query_cluster_index.nil?) and !(hit_cluster_index.nil?)
      # if both hit, merge those 2 clusters into 1
      unless query_cluster_index == hit_cluster_index #do nothing if they are already the same
        # Remove the hit's cluster by moving them all to the query
        clusters[hit_cluster_index].each do |n|
          clusters[query_cluster_index].push n
        end
        clusters.delete_at(hit_cluster_index)
      end
      
    elsif !(query_cluster_index.nil?)
      # if one hits, put the other into the cluster
      clusters[query_cluster_index].push hit
    elsif !(hit_cluster_index.nil?)
      # if one hits, put the other into the cluster
      clusters[hit_cluster_index].push query
    else
      # if neither hits, make a new cluster with the query and hit in it
      clusters.push [query, hit]
    end
    
    
    # record the hit - make a new protein if one doesn't exist, or add a new domain if one already does
    unless options['q']
      dp = domains[hit]
      if dp.nil?
        dp = FramedProtein.new(hit)
        domains[hit] = dp
      end
      dp.add_domain(hit, query_start, query_stop, hit_start, hit_stop)
    else
      dp = domains[query]
      if dp.nil?
        dp = FramedProtein.new(query)
        domains[query] = dp
      end
      dp.add_domain(hit, query_start, query_stop, hit_start, hit_stop)      
    end
    
  end
  
  # Write a file with all the clusters detailed
  clusters.each do |c|
    print c.select{|p| domains[p].nil?}.join(" ")
    print "\t"
    puts c.reject{|p| domains[p].nil?}.join(" ")
  end
  
  # For each cluster, write a chopping file
  clusters.each_with_index do |cluster, i|
    # create a new file in the current directory called clusterX.chopping.csv
    FasterCSV.open("cluster#{i}.chopping.txt",'w',:col_sep => "\t") do |csv|
      # Write out all the extended_hit_domains from the hits to this file
      cluster.each do |name|
        pro = domains[name]
        if pro and !options['q'] #if found a hit and not printing queries
          pro.extended_hit_domains.each do |d|
            csv << [name, d.from2, d.to2]
          end
        elsif pro and options['q'] # if found a query and printing queries
          pro.extended_query_domains.each do |d|
            csv << [name, d.from1, d.to1]
          end
        end
      end
    end
  end
end