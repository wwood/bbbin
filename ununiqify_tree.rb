#!/usr/bin/env ruby

# Take a fasta file and a phylip file that has been uniqified (using for instance
# uniqify_phylip.rb). The only constraint is that they must be in the same order
# in the fasta and the phylip files. Then you rename the nodes of the tree
# (which are currently in phylip unintelligible format) to fasta names,
# which are much more understandable.

require 'bio'

if ARGV.length != 3
  $stderr.puts "Usage: ununiqify_tree.rb <fasta> <uniqued_phylip_file> <tree>"
end

# read the fasta and the phylip files, making a hash between them
fasta_seqs = Bio::FlatFile.open(ARGV[0]).entries
phylip_seqs = Bio::FlatFile.open(Bio::Phylip::PhylipFormat, ARGV[1]).entries[0].alignment.to_fastaformat_array

if fasta_seqs.length != phylip_seqs.size
  $stderr.puts "Number of sequences in fasta and phylip files differ. Are you doing something wrong?"
end

phylip_to_fasta_name_hash = {}
fasta_seqs.each_with_index do |fasta_name, i|
  phylip_to_fasta_name_hash[phylip_seqs[i].definition] = fasta_name.definition
end


# for each node of the tree, rename. warn if there is no hash match
tree = Bio::FlatFile.open(Bio::Newick, ARGV[2]).entries[0].tree
tree.each_node do |node|
  #I get internal nodes here - not sure how else to skip
  next if node.name.nil? or node.name.length == 0 
  
  newname = phylip_to_fasta_name_hash[node.name]
  newname = phylip_to_fasta_name_hash[node.name.gsub(' ','_')] if newname.nil? #bit of a hack
  
  if newname
    node.name = newname
  else
    $stderr.puts "Unexpected node name (left unchanged): '#{node}' '#{node.name}' #{node.class}"
  end
end
puts tree.output(:newick)