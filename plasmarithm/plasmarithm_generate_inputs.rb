#!/usr/bin/env ruby

require 'bio'
require 'bio-plasmoap'
require 'bio-signalp'
require 'reubypathdb'
require 'bio-tm_hmm'
require 'bio-exportpred'
require 'bio-wolf_psort_wrapper'

require 'pp'

# Generated using plasmit_screenscraper.rb
PLASMIT_PREDICTIONS_FILE = '/home/ben/phd/data/Plasmodium falciparum/plasmit/PlasmoDB8.2.plasmit_screenscraped.csv'
# Generated using orthomcl_species_jumper.rb
CRYPTO_ORTHOLOGUES_FILE = 'crypto_orthologues.csv.failed'
# Downloaded directly PlasmoDB
FALCIPARUM_GENE_INFORMATION_PATH = '/home/ben/phd/data/Plasmodium falciparum/genome/PlasmoDB/8.2/PfalciparumGene_PlasmoDB-8.2.txt'
# Supplementary data to the 2006 microarray paper
DERISI_3D7_MICROARRAY_LIFECYCLE_PATH = '/home/ben/phd/data/Plasmodium falciparum/microarray/DeRisi2006/3D7_overview_S6.csv'
# Oehring and Woodcroft et. al. 2012 (hopefully - not yet published)
CORE_NUCLEAR_PROTEOME_PLASMODB_IDS_PATH = '/home/ben/phd/voss_proteome/June2010/nuclear_bioinfo_final6.ids'
# Obtained from Till
HP1_POSITIVE_GENES_PATH = '/home/ben/phd/data/Plasmodium falciparum/HP1chip/hp1.txt'
# Obtained by running falciparum_unconserved_regions.rb
FALCIPARUM_UNCONSERVED_REGION_PLASMODB_IDS = 'falciparum_unconserved_region_genes.csv'

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
  
  # Cache the PlasMit predictions
  falciparum_plasmit_predictions = {}
  File.open(PLASMIT_PREDICTIONS_FILE).each_line do |line|
    splits = line.strip.split("\t")
    raise Exception, "Unexpected PlasMit predictions line format:`#{line}''" if splits.length != 2
    gene_id = Bio::EuPathDB::FastaParser.new('Plasmodium_falciparum_3D7',nil).parse_name(splits[0]).gene_id
    raise Exception, "Found the gene ID `#{gene_id}' multiple times in the P. falciparum plasmit predictions file" if falciparum_plasmit_predictions[gene_id]
    classification = nil
    if splits[1] == 'mito (91%)'
      classification = true
    elsif splits[1] == 'non-mito (99%)'
      classification = false
    else
      raise Exception, "Unexpected Plasmarithm output: `#{splits[1]}'"
    end
    falciparum_plasmit_predictions[gene_id] = classification
  end
  $stderr.puts "Cached #{falciparum_plasmit_predictions.length} Cryptosporidium OrthoMCL orthologues"
  
  # Cache Crypto orthologues file
  falciparum_chom_orthologues = {}
  falciparum_cpar_orthologues = {}
  falciparum_cmur_orthologues = {}
  file_counter = 0
  File.open(CRYPTO_ORTHOLOGUES_FILE).each_line do |line|
    splits = line.strip.split("\t")
    raise unless splits.length >= 2 or splits.length > 5
    raise if falciparum_chom_orthologues[splits[0]]
    falciparum_chom_orthologues[splits[0]] = splits[2]
    falciparum_cpar_orthologues[splits[0]] = splits[3]
    falciparum_cmur_orthologues[splits[0]] = splits[4]
    file_counter += 1
  end
  $stderr.puts "Cached #{file_counter} PlasMit predictions"
  
  # Cache chromosome numbers for each gene
  falciparum_chromosomes = {}
  #  EuPathDBGeneInformationTable.new(File.open(FALCIPARUM_GENE_INFORMATION_PATH)).each do |info|
  #    plasmodb = info.get_info('Gene ID')
  #    chromosome = info.get_info('Chromosome')
  #    raise Exception, "found PlasmoDB ID `#{plasmodb}' multiple times in the Gene information file" if falciparum_chromosomes[plasmodb]
  #    falciparum_chromosomes[plasmodb] = chromosome
  #  end
  #  $stderr.puts "Cached #{falciparum_chromosomes.length} chromosome numbers for P. falciparum genes"
  
  # Cache Bozdech/DeRisi 2006 data
  in_header_section = true
  falciparum_lifecycle_data = {} # Hash of plasmodb => array, where the array hold info each probe corresponding to the gene
  File.open(DERISI_3D7_MICROARRAY_LIFECYCLE_PATH).each_line do |line|
    next if in_header_section and !line.match(/^Oligo/) #skip headers
    in_header_section = false
    next if line.match(/^Oligo/) #skip headers
    
    splits = line.strip.split("\t")
    plasmodb = splits[1]
    info = {
    :amplitude => splits[9],
    :phase => splits[8],
    :timepoint22 => splits[32],
    :timepoint23 => splits[33],
    :timepoint47 => splits[56],
    :timepoint49 => splits[58],
    }
    falciparum_lifecycle_data[plasmodb] ||= []
    falciparum_lifecycle_data[plasmodb].push info
  end
  $stderr.puts "Cached lifecycle microarray data for #{falciparum_lifecycle_data.length} genes"
  
  # Cache nuclear proteome
  nuclear_proteome_plasmodb_ids = []
  File.open(CORE_NUCLEAR_PROTEOME_PLASMODB_IDS_PATH).each_line do |line|
    nuclear_proteome_plasmodb_ids.push line.strip
  end
  
  # Cache HP1 associations
  hp1_positive_genes = []
  File.open(HP1_POSITIVE_GENES_PATH).each_line do |line|
    hp1_positive_genes.push line.strip
  end
  
  # Cache unconserved region genes
  unconserved_region_genes = File.open(FALCIPARUM_UNCONSERVED_REGION_PLASMODB_IDS).readlines.collect{|l| l.strip}

  
  
  
  
  
  
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
    output_line.push Bio::PlasmoAP.new.calculate_score(protein_sequence, signalp.signal?, signalp.cleave(protein_sequence)).points
    
    #  ExportPred?(2): Localisation
    output_line.push !Bio::ExportPred::Wrapper.new.calculate(protein_sequence, :no_RLE => true).predicted_kld?
    output_line.push !Bio::ExportPred::Wrapper.new.calculate(protein_sequence, :no_KLD => true).predicted_kld?
    
    #  WoLF_PSORT prediction Plant(16): Localisation
    #  WoLF_PSORT prediction Animal(15): Localisation
    #  WoLF_PSORT prediction Fungi(12): Localisation
    %w(plant animal fungi).each do |lineage|
      output_line.push Bio::PSORT::WoLF_PSORT::Wrapper.new.run(protein_sequence, lineage).highest_predicted_localization
    end
    
    #  Plasmit(2): Localisation
    if falciparum_plasmit_predictions[plasmodb].nil?
      $stderr.puts "Warning: No PlasMit prediction found for `#{plasmodb}'"
    end
    output_line.push falciparum_plasmit_predictions[plasmodb]
    
    #  Number of C. hominis Genes in Official Orthomcl Group(1): Localisation
    output_line.push !falciparum_chom_orthologues[plasmodb].nil?
    output_line.push !falciparum_cpar_orthologues[plasmodb].nil?
    output_line.push !falciparum_cmur_orthologues[plasmodb].nil?
    
    #  Chromosome(14): Localisation
    output_line.push falciparum_chromosomes[plasmodb]
    
    #  DeRisi 2006 3D7 Timepoint 22(2): Localisation
    #  DeRisi 2006 3D7 Timepoint 23(2): Localisation
    #  DeRisi 2006 3D7 Timepoint 47(2): Localisation
    #  DeRisi 2006 3D7 Timepoint 49(2): Localisation
    if falciparum_lifecycle_data[plasmodb]
      lifecycle_sample = falciparum_lifecycle_data[plasmodb].sample
      [:timepoint22,
       :timepoint23,
       :timepoint47,
       :timepoint49].each do |timepoint|
        output_line.push lifecycle_sample[timepoint]
      end
    else
      output_line.push [nil,nil,nil,nil]
    end
    
    # proteomes
    output_line.push nuclear_proteome_plasmodb_ids.include?(plasmodb)
    
    # number of acidic / basic residue in the first 25 amino acids
    acidics = 0
    protein_sequence[0..24].each_char do |aa|
      acidics += 1 if %w(D E).include?(aa)
    end
    output_line.push acidics
    basics = 0
    protein_sequence[0..24].each_char do |aa|
      basics += 1 if %w(R K H).include?(aa)
    end
    output_line.push basics
    
    # number of transmembrane domains, after the Signal peptide has been removed
    tmhmm = Bio::TMHMM::TmHmmWrapper.new
    result = tmhmm.calculate(signalp.cleave(protein_sequence))
    output_line.push result.transmembrane_domains.length > 0

    # Contained in HP1 associated genomic regions?
    output_line.push hp1_positive_genes.include?(plasmodb)
    
    # In un-aligned regions of the genome (sub-telomeric)?
    output_line.push unconserved_region_genes.include?(plasmodb)

    # dcnls
    # conserved 5' end when blasted against toxo orthologue?
    
    # The answer. What is the localisation?
    
    puts output_line.join(",")
  end
end