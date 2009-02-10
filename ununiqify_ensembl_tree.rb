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
ENSEMBL_SPECIES_HASH = {'ENSRNO'=>['Rat','Rattus norvegicus'],
  'ENSMUS'=>['Mouse','Mus musculus'],
  'ENSGAL'=>['Chicken','Gallus gallus'],
  'ENSBTA'=>['Cow','Bos taurus'],
  'ENSDAR'=>['Zebrafish','Danio rerio'],
  'ENSCAF'=>['Dog','Canis familiaris'],
  'ENSPTR'=>['Chimpanzee','Pan troglodytes'],
  'ENSCPO'=>['Guinea pig','Cavia porcellus'],
  'ENSCIN'=>['C. intestinalis','Ciona intestinalis'],
  'ENSCSAV'=>['C. savignyi','Ciona savignyi'],
  'ENSDNO'=>['Armadillo','Dasypus novemcinctus'],
  'ENSETE'=>['Lesser hedgehog tenrec','Echinops telfairi'],
  'ENSEEU'=>['Hedgehog','Erinaceus europaeus'],
  'ENSFCA'=>['Cat','Felis catus'],
  'ENSGAC'=>['Stickleback','Gasterosteus aculeatus'],
  'ENSLAF'=>['Elephant','Loxodonta africana'],
  'ENSMMU'=>['Macaque','Macaca mulatta'],
  'ENSMOD'=>['Opossum','Monodelphis domestica'],
  'ENSMLU'=>['Microbat','Myotis lucifugus'],
  'ENSOAN'=>['Platypus','Ornithorhynchus anatinus'],
  'ENSOCU'=>['Rabbit','Oryctolagus cuniculus'],
  'ENSORL'=>['Medaka','Oryzias latipes'],
  'ENSSAR'=>['Shrew','Sorex araneus'],
  'ENSSTO'=>['Squirrel','Spermophilus tridecemlineatus'],
  'ENSTBE'=>['Tree shrew','Tupaia belangeri'],
  'SINFRU'=>['Fugu','Takifugu rubripes'],
  'ENSXET'=>['Frog','Xenopus tropicalis'],
  # Below added manually when
  'ENSPCA'=>['Hyrax','Procavia capensis'],
  'ENSOPR'=>['Pika','Ochotona princeps'],
  'ENSDOR'=>['Kangaroo rat','Dipodomys ordii'],
  'ENSECA'=>['Horse','Equus caballus'],
  'ENSPVA'=>['Megabat','Pteropus vampyrus'],
  'ENSTTR'=>['Dolphin','Tursiops truncatus'],
  'ENSMIC'=>['Mouse lemur','Microcebus murinus'],
  'ENSOGA'=>['Bushbaby','Otolemur garnettii'],
  'ENSTRU'=>['Fugu','Takifugu rubripes'],
  'ENSTNI'=>['Tetraodon','Tetraodon nigroviridis'],
  'ENSTSY'=>['Tarsier','Tarsius syrichta'],
  'ENSVPA'=>['Alpaca','Vicugna pacos'],
  'ENSPPY'=>['Orangutan','Pongo pygmaeus'],
  'ENSGGO'=>['Gorilla','gorilla gorilla'],
  'ENS'=>['Human','Homo sapiens']
}
ENSEMBL_OTHER_HASH = {
  /^YEL\d/ => ['Yeast','Saccharomyces cerevisiae'],
  /^T\d/ => ['Worm','Caenorhabditis elegans'],
  /^FBpp\d/ => ['Fly','Drosophila melanogaster'],
  /^AAEL\d/ => ['Aedes','Aedes aegypti'],
  /^AGAP\d/ => ['Anopholes','Anopholes gambiae']
}

class TipLabel
  def initialize(ensembl_name)
    @ensembl_name = ensembl_name
  end
  
  def to_s(use_common_names_only = false)
    # Add in the ensembl name at the beginning
    
    # Try the normal species first
    ENSEMBL_SPECIES_HASH.each do |short, long|
      if @ensembl_name.match(/^#{short}P\d/)
        if use_common_names_only
          return "#{long[0]} #{@ensembl_name.split(/[ \/]/)[0]}"
        else
          return "#{long[0]} (#{long[1]}) #{@ensembl_name.split(/[ \/]/)[0]}"
        end
      end
    end
    
    # If that fails (and therefore doesn't return)
    # Try the other less standard codes
    ENSEMBL_OTHER_HASH.each do |short, long|
      if @ensembl_name.match(short)
        if use_common_names_only
          return "#{long[0]} #{@ensembl_name.split(/[ \/]/)[0]}"
        else
          return "#{long[0]} (#{long[1]}) #{@ensembl_name.split(/[ \/]/)[0]}"
        end
      end 
    end
  
    # Advise if name didn't change
    
    return nil
  end
end

# Below is probably not as DRY as it could be, given that ununiqify_tree.rb does
# something pretty similar. Oh well.

require 'optparse'

# ununiqify_ensembl_tree.rb ../ensembl.cbm48.fa uniqued.phylip consense.outtree
USAGE = "Usage: trees [-c] <fasta_multiple_sequence_alignment> <uniqued.phylip> <consense.outtree>"
options = {
  :common_names => false
}
OptionParser.new do |opts|
  opts.banner = USAGE

  opts.on("-c", "--common-names", "Print common names only, and not scientific names") do |v|
    options[:common_names] = v
  end
end.parse!

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
  
  # set the new name
  newname = TipLabel.new(newname).to_s(options[:common_names])
  
  if newname
    node.name = newname
  else
    $stderr.puts "Unable to find species name for entry id #{newname}"
  end
end
  
puts tree.output(:newick)
