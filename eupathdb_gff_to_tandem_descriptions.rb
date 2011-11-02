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

# Define a class to iterate through the genes in order
class EupathDBGFF
  # Go through each of the genes in the GFF file, after
  # ordering them
  def ordered_gene_iterator(&block)
    # First cache all the genes
    genes = []
    while (g = next_gene)
      genes.push g
    end
    
    # Then sort them according to chromosome and start point
    genes.sort! do |g1, g2|
      chr_match = (g1.seqname <=> g2.seqname)
      if chr_match != 0
        chr_match
      else
        g1.start <=> g2.start #compare starts. Could be more specific here but there are probably no cases where 2 genes have the same start
      end      
    end
    
    # Then yield in order
    genes.each do |g|
      yield g
    end
  end
end

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

chromosomes_annotated = {}

# Don't analyse mitochondrial or apicoplast genomes, for instance
uninterrogated_scaffolds = [
  'apidb|PFC10_API_IRAB'
]

gff.ordered_gene_iterator do |current_gene|
  raise if current_gene.nil? #should never happen if the iterator is working
  num += 1
  
  # ignore when they cross chromosome boundaries, and when they are on non-nuclear chromosomes
  if current_gene.seqname == last_gene.seqname and !uninterrogated_scaffolds.include?(current_gene.seqname)
   
  #debug 
  #$stderr.puts [last_gene.name, current_gene.name, last_gene.strand, current_gene.strand].join("\t")
  
  # Go through the OrthoMCL list
  current_orthomcl = get_orthomcl_group.call current_gene
  
  # What is the butting?
  butting = nil
  if last_gene.strand == '+'
    if current_gene.strand == '+'
      butting = 'tail_to_head'
    elsif current_gene.strand == '-'
      butting = 'tail_to_tail'
    else
      raise
    end
  elsif last_gene.strand == '-'
    if current_gene.strand == '+'
      butting = 'head_to_head'
    elsif current_gene.strand == '-'
      butting = 'head_to_tail'
    else
      raise
    end
  else
    raise
  end
  
  chromosomes_annotated[current_gene.seqname] ||= 0
  chromosomes_annotated[current_gene.seqname] += 1
  
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
  #break if num > 100
end

chromosomes_annotated.each do |name, num|
  $stderr.puts [name, num].join("\t")
end
