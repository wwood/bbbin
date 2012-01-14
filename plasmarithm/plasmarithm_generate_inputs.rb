#!/usr/bin/env ruby

require 'bio'
require 'bio-plasmoap'
require 'bio-signalp'
require 'reubypathdb'
require 'bio-tm_hmm'
require 'bio-exportpred'
require 'bio-wolf_psort_wrapper'

if __FILE__ == $0
  # Take a file that contains 1 PlasmoDB ID per line, and generate a matrix that 
  # corresponds to the inputs for use in R.
  
  USAGE = 'plasmarithm_generate_inputs.rb <plasmodb_id_list_file>'
  
  # First, cache the protein sequences
  falciparum = EuPathDBSpeciesData.new('Plasmodium falciparum','/home/ben/phd/data')
  falciparum_protein_sequences = {}
  falciparum.protein_fasta_file_iterator.each do |prot|
    gene_id = prot.gene_id
    raise Exception, "Found the gene ID `#{gene_id}' multiple times in the P. falciparum protein fasta file" if falciparum_protein_sequences[gene_id]
    falciparum_protein_sequences[gene_id] = prot.sequence
  end
  $stderr.puts "Cached #{falciparum_protein_sequences.length} protein sequences"
  
  
  # For each PlasmoDB ID provided in the input file
  ARGF.each do |line|
    plasmodb = line.strip
    protein_sequence = falciparum_protein_sequences[plasmodb]
    if protein_sequence.nil?
      raise Exception, "Unable able to find protein sequence for PlasmoDB ID `#{plasmodb}'"
    end
    
    # This gets progressively filled with data about the current gene
    output_line = []
    
    #  SignalP Prediction(2): Localisation
    signalp = Bio::SignalP::Wrapper.new.calculate(protein_sequence)
    output_line.push signalp.signal?
    
    #  PlasmoAP Score(2): Localisation
    # We already know whether there is a SignalP prediction or not, so just go with that. Re-calculating takes time.
    plasmoap = Bio::PlasmoAP.new.calculate_score(protein_sequence, signalp.signal?, signalp.cleave(protein_sequence))
    
    #  ExportPred?(2): Localisation
    output_line.push Bio::ExportPred::Wrapper.new.calculate(protein_sequence, :no_KLD => true).signal?
    output_line.push Bio::ExportPred::Wrapper.new.calculate(protein_sequence, :no_RLE => true).signal?
    
    #  WoLF_PSORT prediction Plant(16): Localisation
    #  WoLF_PSORT prediction Animal(15): Localisation
    #  WoLF_PSORT prediction Fungi(12): Localisation
    %w(plant animal fungi).each do |lineage|
      output_line.push Bio::PSORT::WoLF_PSORT.new.run(protein_sequence, lineage).highest_predicted_localization
    end

    #  Plasmit(2): Localisation
    

    #  Number of C. homonis Genes in Official Orthomcl Group(1): Localisation
    #  Chromosome(14): Localisation
    #  DeRisi 2006 3D7 Timepoint 22(2): Localisation
    #  DeRisi 2006 3D7 Timepoint 23(2): Localisation
    #  DeRisi 2006 3D7 Timepoint 47(2): Localisation
    #  DeRisi 2006 3D7 Timepoint 49(2): Localisation 
    
    # proteomes
    # number of acidic / basic residue in the first 25 amino acids
    # number of transmembrane domains
    # Contained in HP1 associated genomic regions?
    
    puts output_line.join(",")
  end
end