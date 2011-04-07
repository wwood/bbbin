#!/usr/bin/env ruby

# A script to make the design primers for T. gondii ku80 3' replacement. Probably
# doesn't work until further notice.

require 'rubygems'
require 'bio'
require 'optparse'
require 'eupathdb_gene_information_table' #from the reubypathdb rubygem



if $0 == __FILE__
  options = {
    :gene_information_table_path => '/home/ben/phd/data/Toxoplasma gondii/ToxoDB/6.3/TgondiiME49Gene_ToxoDB-6.3.txt'
  }
  wanted_gene_id = nil
  o = OptionParser.new do |opts|
    opts.banner = [
      'Usage: program <toxodb_ME49_id>'
    ]
    if ARGV.length != 1
    	$stderr.puts opts.banner
    	exit
    end
  end
o.parse!

wanted_gene_id = ARGV[0]



# Delete this later
# A class for extracting gene info from a gene info file
class EuPathDBGeneInformationFileExtractor
  def initialize(filename = nil)
    @filename = filename
  end

  # Returns a EuPathDBGeneInformation object corresponding to the wanted key. If
  # there are multiple in the file, only the first is returned.
  def extract_gene_info(wanted_gene_id)
    EuPathDBGeneInformationTable.new(File.open(@filename)).each do |gene|
      return gene if wanted_gene_id == gene.get_info('Gene Id')
    end
    return nil
  end
end


# Input the position that is to be inserted. Presumably a ToxoDB Gene id, assuming a 3' insertion at the end of the gene just before the stop codon
# First, extract gene info
gene_info = EuPathDBGeneInformationFileExtractor.new(options[:gene_information_table_path]).extract_gene_info(wanted_gene_id)

# Parse +/-ve direction of gene
location = gene_info.get_info('Genomic Location')
$stderr.puts "Genomic location: #{location}"
direction = nil #true => fwd, false => backwards
chromosome = nil
if matches = location.match(/(.*): [\d,]+ - [\d,]+ \(([+-])\)$/)
	chromosome = matches[1]
	direction = matches[2]
else
	raise Exception, "Couldn't parse genomic location line `#{location}'"
end
$stderr.puts "Found a #{direction} direction sequence on #{chromosome}"

# Parse out the end exon - highest ending if +ve direction, lowest ending if -ve direction
# It doesn't appear to be possible to extract the cds positions from the info file, only the transcript exons'

# Extract 1kb upstream of insertion point



# Find the unique restriction sites that are in that 1kb
# Remove the restriction sites that are known to be in the plasmid we are inserting into
# Choose one of the enzymes to cut with. It must be at least 350bp from the 3' end of the sequence

# Attempt to design primers
# right primer: using the reverse complement of the first however many bases
# left primer: try to design this
# ensure that the length of the product is at least the distance from the  end of the sequence to the restriction site + 350bp

# If primers found, output, otherwise
# a) try another restriction enzyme
# b) try with more than 1kb being extracted

end