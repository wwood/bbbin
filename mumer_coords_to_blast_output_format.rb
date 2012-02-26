#!/usr/bin/env ruby

#Convert a MUMER coords file output to something like blast

# /srv/whitlam/home/users/uqbwoodc/ABISKO/Split_Library_Output2/uclust_picked_otus/rep_set.fna /srv/whitlam/home/users/uqbwoodc/metagenome2/mira/assemblies/1/analyses/rep_set_then_greengenes/../metagenome2_out.unpadded.fasta
# NUCMER
#
# [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  | [TAGS]
# =====================================================================================
# 1      186  |      199      384  |      186      186  |   100.00  | 164  metagenome2_c17478

indices = {}

ARGF.each_line do |line|
  if $. == 4
    headers = line.strip.split(/[ \|]\s+/)
    {
      '[S1]' => :query_start,
      '[S2]' => :query_end,
      '[TAGS]' => :query_name
    }.each do |header, sym|
      headers.each_with_index do |head, i|
        if head == header
          indices[sym] = i
        end
      end
      raise Exception, "Unable to find index for #{sym}" if indices[sym].nil?
    end
  end
  
  next unless $. > 5 #skip header bits
  splits = line.strip.split(/\s+/)

  puts [
    splits[indices[:query_name]],
    %w(. . . . . ),
    splits[indices[:query_start]],
    splits[indices[:query_end]],
  ].join("\t")
end
