#!/usr/bin/env ruby

#Convert a MUMER coords file output to something like blast

# /srv/whitlam/home/users/uqbwoodc/ABISKO/Split_Library_Output2/uclust_picked_otus/rep_set.fna /srv/whitlam/home/users/uqbwoodc/metagenome2/mira/assemblies/1/analyses/rep_set_then_greengenes/../metagenome2_out.unpadded.fasta
# NUCMER
#
# [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  | [TAGS]
# =====================================================================================
# 1      186  |      199      384  |      186      186  |   100.00  | 164  metagenome2_c17478

ARGF.each_line do |line|
  next unless $. > 5 #skip header bits
  splits = line.strip.split(/\s+/)
  raise unless splits.length == 19
  puts [
    splits[17],
    %w(. . . . . ),
    splits[0],
    splits[3],
  ].join("\t")
end
