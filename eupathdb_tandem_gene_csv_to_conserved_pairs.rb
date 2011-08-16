#!/usr/bin/env ruby

require 'pp'

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
  next if o1.nil? or o2.nil?
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
  next if o1.nil? or o2.nil?
  key = [o1,o2].sort
  
  orthomcl_to_pairs2[key] ||= []
  orthomcl_to_pairs2[key].push [gene1, gene2, butting]
end

# pp orthomcl_to_pairs2.length
# pp orthomcl_to_pairs2.to_a[0]
# pp orthomcl_to_pairs2.to_a[1]


# Find conserved gene pairs
orthomcl_to_pairs1.each do |one_orths, ones|
  next unless orthomcl_to_pairs1[one_orths].length == 1 #only count 1 to 1 orthologues
  twos_array = orthomcl_to_pairs2[one_orths]
  if !twos_array.nil? and twos_array.length == 1
    twos = twos_array[0]
    puts [one_orths,ones,twos].flatten.join("\t")
  end
end