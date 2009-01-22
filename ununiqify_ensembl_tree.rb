#!/usr/bin/env ruby

# Take a fasta file and a phylip file that has been uniqified (using for instance
# uniqify_phylip.rb). The only constraint is that they must be in the same order
# in the fasta and the phylip files. Then you rename the nodes of the tree
# (which are currently in phylip unintelligible format) to fasta names,
# which are much more understandable.
#
# Assumes also that the names in the fasta file are named at least in the beginning
# like Ensembl genes, and then prepends the species name to the front, so that
# things are more obvious in a tree

require 'bio'

# Ensembl species taken from ensembl/modules/Bio/EnsEMBL/Registry.pm
# Ensembl version 52 (December 2009)
# transformed manually to hash form
ENSEMBL_SPECIES_HASH = {'ENSRNO'=>'Rattus norvegicus',
  'ENSMUS'=>'Mus musculus',
  'ENSGAL'=>'Gallus gallus',
  'ENSBTA'=>'Bos taurus',
  'ENSDAR'=>'Danio rerio',
  'ENSCAF'=>'Canis familiaris',
  'ENSPTR'=>'Pan troglodytes',
  'ENSCPO'=>'Cavia porcellus',
  'ENSCIN'=>'Ciona intestinalis',
  'ENSCSAV'=>'Ciona savignyi',
  'ENSDNO'=>'Dasypus novemcinctus',
  'ENSETE'=>'Echinops telfairi',
  'ENSEEU'=>'Erinaceus europaeus',
  'ENSFCA'=>'Felis catus',
  'ENSGAC'=>'Gasterosteus aculeatus',
  'ENSLAF'=>'Loxodonta africana',
  'ENSMMU'=>'Macaca mulatta',
  'ENSMOD'=>'Monodelphis domestica',
  'ENSMLU'=>'Myotis lucifugus',
  'ENSOAN'=>'Ornithorhynchus anatinus',
  'ENSOCU'=>'Oryctolagus cuniculus',
  'ENSORL'=>'Oryzias latipes',
  'ENSSAR'=>'Otolemur garnettii',
  'ENSSTO'=>'Spermophilus tridecemlineatus',
  'ENSTBE'=>'Tupaia belangeri',
  'SINFRU'=>'Takifugu rubripes',
  'ENSXET'=>'Xenopus tropicalis'
}
ENSEMBLE_HUMAN_HASH = {'ENS'=>'Homo sapiens'}


# Below is probably not as DRY as it could be, given that ununiqify_tree.rb does
# something pretty similar. Oh well.

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
    next
  end
  
  # Add in the ensembl name at the beginning
  changed = false
  ENSEMBL_SPECIES_HASH.each do |short, long|
    if newname.match(/^#{short}/)
      changed = true
      newname = "#{long} #{newname.split(' ')[0]}"
    end
  end
  unless changed # Only do human if I have to
    ENSEMBLE_HUMAN_HASH.each do |short, long|
      if newname.match(/^#{short}/)
        changed = true
        newname = "#{long} #{newname.split(' ')[0]}"
      end
    end
  end
  
  # Advise if name didn't change
  if changed
    
  else
    newname = newname.split(' ')[0]
    $stderr.puts "Unable to find species name for entry id #{newname}"
  end
  node.name = newname
end
  
puts tree.output(:newick)
