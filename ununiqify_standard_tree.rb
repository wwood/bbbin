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

require 'rubygems'
require 'bio'
require File.dirname(__FILE__) + '/ensembl'

class TipLabel

  def initialize(ensembl_name)
    @ensembl_name = ensembl_name
  end
  
  def to_s(use_common_names_only = false)
    # Add in the ensembl name at the beginning
    
    # Try the normal species first
    Bio::Ensembl::ENSEMBL_SPECIES_HASH.each do |short, long|
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
    Bio::Ensembl::ENSEMBL_OTHER_HASH.each do |short, long|
      if @ensembl_name.match(short)
        if use_common_names_only
          return "#{long[0]} #{@ensembl_name.split(/[ \/]/)[0]}"
        else
          return "#{long[0]} (#{long[1]}) #{@ensembl_name.split(/[ \/]/)[0]}"
        end
      end
    end

    # Try JGI tree
    Bio::JGI::JGI_SPECIES_HASH.each do |short, long|
      if @ensembl_name.match(/^#{short}\d\|/)
        if use_common_names_only
          return "#{long[0]} #{@ensembl_name.split(/[ \/]/)[0]}"
        else
          return "#{long[0]} (#{long[1]}) #{@ensembl_name.split(/[ \/]/)[0]}"
        end
      end
    end

    # Try other species specified as regexs
    [Bio::EuPathDB::EUPATHDB_SPECIES_HASH, Bio::Misc::MISC_SPECIES_HASH].each do |hash|
      hash.each do |short, long|
        if @ensembl_name.match(short)
          if use_common_names_only
            return "#{long[0]} #{@ensembl_name.split(/[ \/]/)[0]}"
          else
            return "#{long[0]} (#{long[1]}) #{@ensembl_name.split(/[ \/]/)[0]}"
          end
        end
      end
    end

    # Try Orthomcl, which doesn't have common names
    [Bio::Orthomcl::ORTHOMCL_SPECIES_HASH].each do |hash|
      hash.each do |short, long|
        if @ensembl_name.match(short)
          if use_common_names_only
            return "#{long[0]} #{@ensembl_name.split(/[ \/]/)[0]}"
          else
            return "#{long[0]} #{@ensembl_name.split(/[ \/]/)[0]}"
          end
        end
      end
    end
  
    # Advise if name didn't change
    $stderr.puts "Name not recognized - left unchanged: '#{@ensembl_name}'"
    return @ensembl_name
  end
end

# Below is probably not as DRY as it could be, given that ununiqify_tree.rb does
# something pretty similar. Oh well.

require 'optparse'

# ununiqify_ensembl_tree.rb ../ensembl.cbm48.fa uniqued.phylip consense.outtree
USAGE = "Usage: ununiqify_ensembl_tree [-c] [-m <manual_names_filename>] <fasta_multiple_sequence_alignment> <uniqued.phylip> <consense.outtree>"
options = {
  :common_names => false,
  :manual_names => {},
  :phylip_manual_names = {},
}
OptionParser.new do |opts|
  opts.banner = USAGE

  opts.on("-c", "--common-names", "Print common names only, and not scientific names") do |v|
    options[:common_names] = v
  end

  opts.on("-m", "--manual-names MANUAL_NAMES_FILENAME", "Read a hash of regular expression => wanted names for a node, so that figures can be recreated more easily (may be BUGGGGY)") do |v|
    require v
    options[:manual_names] = MANUAL_NAMES
  end

  opts.on("-p", "--phylip-manual-names MANUAL_NAMES_FILENAME", "Read a hash of names in the uniqued phylip file => wanted names for a node, so that figures can be recreated more easily") do |v|
    require v
    options[:phylip_manual_names] = PHYLIP_MANUAL_NAMES
  end
end.parse!

unless ARGV.length == 3
  $stderr.puts USAGE
  exit 1
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
phylip_to_fasta_name_hash = {}

# for each node of the tree, rename. warn if there is no hash match
tree = Bio::FlatFile.open(Bio::Newick, ARGV[2]).entries[0].tree
tree.leaves.each do |node|

  # Priority of renaming:
  # 0. manually specified names of the uniqued phylip file
  # 1. manually specified as a string in the manual_names hash, where the key is the fasta filename
  # 2. manually specified as a regular expression in the manual_names hash
  # 3. using TipLabel to rename common species
  # 4. Use the original name (and warn that this is happening)

  # 3. happens first because this is what is used in the manual matching
  newname = phylip_to_fasta_name_hash[node.name]
  newname = phylip_to_fasta_name_hash[node.name.gsub(' ','_')] if newname.nil? #bit of a hack

  manualled = false

  # 0. manually specified names of the uniqued phylip file
  options[:phylip_manual_names].each do |key, value|
    next unless key.kind_of?(String)
    if key == node.name
      newname = value
      manualled = true
    end
  end

  # 1. manually specified as a string in the manual_names hash
unless manualled
  options[:manual_names].each do |key, value|
    next unless key.kind_of?(String)
    if key == newname
      newname = value
      manualled = true
    end
  end
end

  # 2. manually specified as a regular expression in the manual_names hash
  unless manualled
    options[:manual_names].each do |key, value|
      next unless key.kind_of?(Regexp)
#      $stderr.puts key.inspect
#      $stderr.puts newname
      if newname.match(key)
        newname = value
        manualled = true
      end
    end
  end

  original = newname
  
  if newname
    node.name = newname
  else
    # 4. Use the original name (and warn that this is happening)
    $stderr.puts phylip_to_fasta_name_hash.inspect
    $stderr.puts "Unexpected node name (left unchanged): '#{node}' '#{node.name}' #{node.class}"
    next
  end
  
  # set the new name in the tree
  if manualled
    node.name = newname
  elsif newname
    node.name = TipLabel.new(newname).to_s(options[:common_names])
  else
    $stderr.puts "Unable to find species name for entry id '#{original}'"
  end
end
  
puts tree.output(:newick)
