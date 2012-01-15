#!/usr/bin/env ruby


# Find groups of genes in the P. falciparum genome that are not conserved to other Plasmodium species


require 'reubypathdb'
require 'pp'

if __FILE__ == $0
  FALCIPARUM_GENE_INFORMATION_FILE = '/home/ben/phd/data/Plasmodium falciparum/genome/PlasmoDB/8.2/PfalciparumGene_PlasmoDB-8.2.txt'
  FALCIPARUM_ONE_TO_ONE_ORTHOLOGUES_FILE = '/home/ben/phd/amino_acid_overrepresentation/1/plasmodium_orthologues/pfal_pviv_pber_pcha.pfal.genes'
  
  WINDOW_SIZE=10
  MAX_ONE_TO_ONE_ORTHOLOGUES_PER_WINDOW=1
  
  falciparum_chromosomes = {} #Hash of chromosome => plasmodb_id => gene start 
  count = 0
  EuPathDBGeneInformationTable.new(File.open(FALCIPARUM_GENE_INFORMATION_FILE)).each do |info|
    plasmodb = info.get_info('Gene ID')
    chromosome = info.get_info('Chromosome')
    next if info.get_info('CDS Length') == 'null'
    next unless (1..14).collect{|i| i.to_s}.include?(chromosome)
    location = info.get_info('Genomic Location')
    # Genomic Location: Pf3D7_14: 778,091 - 778,162 (+)
    if matches = location.match(/: ([\d,]+) - ([\d,]+)/)
      one = matches[1].gsub(',','').to_i
      two = matches[2].gsub(',','').to_i
      raise if one == 0 or two == 0
      first = one<two ? one : two
      
      falciparum_chromosomes[chromosome] ||= {}
      falciparum_chromosomes[chromosome][plasmodb] = first
    else
      $stderr.puts "Couldn't parse genomic location for #{plasmodb}: #{location}"
    end
    count += 1
    #break if count > 20 #debug
  end
  $stderr.puts "Cached #{count} P. falciparum genes from the gene info file"
  
  # Cache which genes have 1:1 orthology
  one_to_one_orthology = File.open(FALCIPARUM_ONE_TO_ONE_ORTHOLOGUES_FILE).readlines.collect{|l| l.strip}
  
  falciparum_chromosomes.each do |chromosome, genes|
    current_window = []
    unconserved = []
    genes.to_a.sort{|a,b|
      a[1] <=> b[1]
    }.each do |gene|
      # If we have a full window here
      if current_window.length == WINDOW_SIZE
        # First, is the current window unconserved?
        num_conserved = current_window.select {|g|
          one_to_one_orthology.include?(g[0])
        }.length
        if num_conserved <= MAX_ONE_TO_ONE_ORTHOLOGUES_PER_WINDOW
          current_window.each{|g| unconserved.push g[0] unless one_to_one_orthology.include?(g[0])}
        end
        current_window = current_window[1..(WINDOW_SIZE-1)]
      end
      current_window.push gene
    end
    puts unconserved.uniq.join("\n")
  end
end