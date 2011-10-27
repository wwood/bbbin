#!/usr/bin/env ruby

require 'pp'
# Available from http://mbe.oxfordjournals.org/content/suppl/2011/04/19/msr103.DC1/DeBarry_SuppTable1_MBE-10-1053.txt
debarry_supp_path = "/home/ben/phd/data/apicomplexa/synteny/DeBarry_SuppTable1_MBE-10-1053.txt"

csv1 = ARGV[0]
csv2 = ARGV[1]

orthomcl_to_pairs1 = {}
orthomcl_to_pairs2 = {}

File.open(csv1).each_line do |a|
  things = a.chomp.split("\t")
  gene1 = things[0]
  gene2 = things[1]
  butting = things[2]
  o1 = things[3]
  o2 = things[4]
  next if o1.nil? or o2.nil? or o1=="" or o2==""
  key = [o1,o2].sort
  
  orthomcl_to_pairs1[key] ||= []
  orthomcl_to_pairs1[key].push [gene1, gene2, butting]
end

# pp orthomcl_to_pairs1.length
# pp orthomcl_to_pairs1.to_a[0]
# pp orthomcl_to_pairs1.to_a[1]

File.open(csv2).each_line do |a|
  things = a.chomp.split("\t")
  gene1 = things[0]
  gene2 = things[1]
  butting = things[2]
  o1 = things[3]
  o2 = things[4]
  next if o1.nil? or o2.nil? or o1=="" or o2==""
  key = [o1,o2].sort
  
  orthomcl_to_pairs2[key] ||= []
  orthomcl_to_pairs2[key].push [gene1, gene2, butting]
end

# pp orthomcl_to_pairs2.length
# pp orthomcl_to_pairs2.to_a[0]
# pp orthomcl_to_pairs2.to_a[1]

# Assign each gene as a syntenic block or not
syntenic_gene_pairs = {} #hash of gene pairs to block number
syntenic_genes = []
File.open(debarry_supp_path).each_line do |line|
  next if line.match /^#/
  line.chomp!
  splits = line.split ' '
  next unless splits.length == 9
  syntenic_gene_pairs[[splits[2],splits[5]].sort] ||= []
  syntenic_gene_pairs[[splits[2],splits[5]].sort].push splits[0]
end

syntenic_gene_pairs.to_a[0..10].each do |a|
  #$stderr.puts a.join(', ')
end

# Find conserved gene pairs
puts [
  'OrthoMCL_1',
  'OrthoMCL_2',
  'species_1_gene_1',
  'species_1_gene_2',
  'species_1_butting',
  'gene_1_synteny',
  'species_2_gene_1',
  'species_2_gene_2',
  'species_2_butting',
  'gene_2_synteny',
].join("\t")
orthomcl_to_pairs1.each do |one_orths, ones|
  next unless orthomcl_to_pairs1[one_orths].length == 1 #only count 1 to 1 orthologues
  twos_array = orthomcl_to_pairs2[one_orths]
  if !twos_array.nil? and twos_array.length == 1
    twos = twos_array[0]
    syntenic1 = syntenic_gene_pairs[[twos[0],ones[0][1]].sort]
    syntenic2 = syntenic_gene_pairs[[ones[0][0],twos[1]].sort]
    syntenic1 ||= []; syntenic2 ||= [];
    puts [one_orths,ones,syntenic1.join(','),twos,syntenic2.join(',')].flatten.join("\t")
  end
end