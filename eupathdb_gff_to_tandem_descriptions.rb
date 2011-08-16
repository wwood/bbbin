#!/usr/bin/env ruby

# Take a EuPathDB format GFF file, and convert it into a file
# describing tandem genes and map each gene to their corresponding OrthoMCL group

# Assumes the genes are in order in the GFF file

require 'rubygems'
require 'reubypathdb'
require 'bio'
require 'bio-orthomcl'
require 'pp'
require 'zlib'

unless ARGV.length == 3
  $stderr.puts "Usage: #{$0} orthomcl_groups_path species_code gff_file_path"
  exit 1
end

orthomcl_groups_path = ARGV[0]
species_code = ARGV[1]
gff_file_path = ARGV[2]

# First, cache the OrthoMCL file's contents as a hash
io = Zlib::GzipReader.open(orthomcl_groups_path)
orthomcl_hash = Bio::OrthoMCL.convert_groups_file_to_gene_id_to_group_id_hash(io, species_code)

# Go through the GFF file, looking for tandem genes
gff = EupathDBGFF.new(gff_file_path)
last_gene = gff.next_gene
raise if last_gene.nil?
next_gene = nil

get_orthomcl_group = lambda do |gene|
  group = nil
  [gene.name, gene.alternate_ids].flatten.each do |i|
    possible = "#{species_code}|#{i}"
    group = orthomcl_hash[possible]
    break unless group.nil?
  end
  group
end
last_orthomcl = get_orthomcl_group.call last_gene
num = 0

while true
  # alternate_ids is the list of alternate names, and these may be used in OrthoMCL
  current_gene = gff.next_gene
  break if current_gene.nil? #escape the infinite loop
  num += 1
  
  # ignore when they cross chromosome boundaries
  if current_gene.seqname == last_gene.seqname
  
  # Go through the OrthoMCL list
  current_orthomcl = get_orthomcl_group.call current_gene
  
  # What is the butting?
  butting = nil
  if last_gene.strand == '+'
    if current_gene.strand = '+'
      butting = 'tail_to_head'
    elsif current_gene.strand = '-'
      butting = 'tail_to_tail'
    else
      raise
    end
  elsif last_gene.strand == '-'
    if current_gene.strand = '+'
      butting = 'head_to_head'
    elsif current_gene.strand = '-'
      butting = 'head_to_tail'
    else
      raise
    end
  else
    raise
  end
  
  puts [
  last_gene.name,
  current_gene.name,
  butting,
  last_orthomcl,
  current_orthomcl,
  ].join("\t")
  end
  
  last_gene = current_gene
  last_orthomcl = current_orthomcl
  #break if num > 10
end
