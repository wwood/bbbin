#!/usr/bin/env ruby

require 'bio'
require 'bio-plasmoap'
require 'bio-signalp'
require 'reubypathdb'
require 'bio-tm_hmm'
require 'bio-exportpred'
require 'bio-wolf_psort_wrapper'

require 'pp'
require 'progressbar'

# Generated using plasmit_screenscraper.rb
PLASMIT_PREDICTIONS_FILE = '/home/ben/phd/data/Plasmodium falciparum/plasmit/PlasmoDB8.2.plasmit_screenscraped.csv'
# Generated using orthomcl_species_jumper.rb
CRYPTO_ORTHOLOGUES_FILE = 'crypto_orthologues.csv.failed'
# Downloaded directly PlasmoDB
FALCIPARUM_GENE_INFORMATION_PATH = '/home/ben/phd/data/Plasmodium falciparum/genome/PlasmoDB/8.2/PfalciparumGene_PlasmoDB-8.2.txt'
# Supplementary data to the 2006 microarray paper
DERISI_3D7_MICROARRAY_LIFECYCLE_PATH = '/home/ben/phd/data/Plasmodium falciparum/microarray/DeRisi2006/3D7_overview_S6.csv'
# Oehring and Woodcroft et. al. 2012 (hopefully - not yet accepted!)
CORE_NUCLEAR_PROTEOME_PLASMODB_IDS_PATH = '/home/ben/phd/voss_proteome/June2010/nuclear_bioinfo_final6.ids'
# Obtained from Till
HP1_POSITIVE_GENES_PATH = '/home/ben/phd/data/Plasmodium falciparum/HP1chip/hp1.txt'
# Obtained by running falciparum_unconserved_regions.rb
FALCIPARUM_UNCONSERVED_REGION_PLASMODB_IDS = 'falciparum_unconserved_region_genes.window5max1.csv'
# used dcnls.rb and some other things for these two
DCNLS_4_OF_5_PATH = 'dcnls.4of5.csv'
DCNLS_5_OF_6_PATH = 'dcnls.5of6.csv'
# Conservation of 5' end of the gene
FALCIPARUM_TO_TOXO_BLAST_CSV = 'falciparum8.2versusToxo8.2.blast.csv'



class Array
    def sum
        self.inject{|sum,x| sum + x }
    end
end



if __FILE__ == $0
  # Take a file that contains 1 PlasmoDB ID per line, and generate a matrix that 
  # corresponds to the inputs for use in the R package glmnet.
  
  USAGE = 'plasmarithm_generate_inputs.rb <plasmodb_id_list_file>'
  
  # Cache the falciparum->toxo prediction
  falciparum_to_toxo_min_starts = {}
  File.open(FALCIPARUM_TO_TOXO_BLAST_CSV).each_line do |line|
    splits = line.strip.split("\t")
    pl = splits[0].gsub(/^..../,'') #get rid of the 'psu|' at the start
    # Fields: 0query id, 1subject id, 2% identity, 3alignment length, 4ismatches, 5gap opens, 6q. start, q. end, s. start, s. end, evalue, bit score
    start = splits[6].to_i
    if falciparum_to_toxo_min_starts[pl].nil? or falciparum_to_toxo_min_starts[pl]>start
      falciparum_to_toxo_min_starts[pl] = start
    end 
  end
  
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
  EuPathDBGeneInformationTable.new(File.open(FALCIPARUM_GENE_INFORMATION_PATH)).each do |info|
    plasmodb = info.get_info('Gene ID')
    chromosome = info.get_info('Chromosome')
    raise Exception, "found PlasmoDB ID `#{plasmodb}' multiple times in the Gene information file" if falciparum_chromosomes[plasmodb]
    falciparum_chromosomes[plasmodb] = chromosome
  end
  $stderr.puts "Cached #{falciparum_chromosomes.length} chromosome numbers for P. falciparum genes"
  
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
    }
    falciparum_lifecycle_data[plasmodb] ||= []
    falciparum_lifecycle_data[plasmodb].push info
    info[:timepoints] = []
    (1..53).each do |i|
      meas = splits[i+9]
      meas = 1 if meas=='NA'
      meas = meas.to_f
      info[:timepoints].push meas
    end
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
  
  # Cache DCNLS output files
  dcnls4of5 = {}
  File.open(DCNLS_4_OF_5_PATH).readlines.each do |l|
    splits = l.split("\t")
    plasmodb = splits[1]
    num = splits[2].length > 0 ? splits[2].split(',') : nil
    dcnls4of5[plasmodb] = num
  end
  dcnls5of6 = {}
  File.open(DCNLS_5_OF_6_PATH).readlines.each do |l|
    splits = l.split("\t")
    plasmodb = splits[1]
    num = splits[2].length > 0 ? splits[2].split(',') : nil
    dcnls5of6[plasmodb] = num
  end
  

  
  
  
  
  
  
  # For each PlasmoDB ID provided in the input file
  do_headers = true
  headers = []
  plasmodbs = ARGF.readlines
  progress = ProgressBar.new('inputs', plasmodbs.length)
  plasmodbs.each do |line|
    plasmodb = line.strip
    protein_sequence = falciparum_protein_sequences[plasmodb]
    if protein_sequence.nil?
      raise Exception, "Unable able to find protein sequence for PlasmoDB ID `#{plasmodb}'"
    end
    
    # This gets progressively filled with data about the current gene
    output_line = []
    
    headers.push 'plasmodb' if do_headers
    output_line.push plasmodb
    
    #  SignalP Prediction(2): Localisation
    headers.push 'signalp' if do_headers
    signalp = Bio::SignalP::Wrapper.new.calculate(protein_sequence)
    output_line.push signalp.signal?
    
    #  PlasmoAP Score(2): Localisation
    # We already know whether there is a SignalP prediction or not, so just go with that. Re-calculating takes time.
    headers.push 'plasmoap' if do_headers
    output_line.push Bio::PlasmoAP.new.calculate_score(protein_sequence, signalp.signal?, signalp.cleave(protein_sequence)).points
    
    #  ExportPred?(2): Localisation
    headers.push ['exportpredRLE', 'exportpredKLD'] if do_headers
    output_line.push !Bio::ExportPred::Wrapper.new.calculate(protein_sequence, :no_RLE => true).predicted_kld?
    output_line.push !Bio::ExportPred::Wrapper.new.calculate(protein_sequence, :no_KLD => true).predicted_kld?
    
    #  WoLF_PSORT prediction Plant(16): Localisation
    #  WoLF_PSORT prediction Animal(15): Localisation
    #  WoLF_PSORT prediction Fungi(12): Localisation
    # Encode each with the scores that they were given, accounting for dual locations
    wolfer = lambda do |lineage, all_locations|
      all_locations.collect{|l| l.downcase}.sort.each{|l| headers.push "wolf_#{lineage}_#{l}"} if do_headers
      result = Bio::PSORT::WoLF_PSORT::Wrapper.new.run(protein_sequence, lineage)
      
      # Initialise 
      location_scores = {}
      all_locations.each {|l| location_scores[l] = 0.0}
      result.score_hash.each do |loc, score|
        score = score.to_f
        raise unless score >0.0
        if matches = loc.match(/(.+)_(.+)/) #if dual localised, add the score to both locations
          raise unless location_scores[matches[1]] and location_scores[matches[2]]
          location_scores[matches[1]] += score/2.0
          location_scores[matches[2]] += score/2.0
        else #regular 1 location prediction
          raise Exception, "WoLF PSORT predicted an unexpected #{lineage} location: `#{loc}'" unless location_scores[loc]
          location_scores[loc] += score
        end
      end
      # Write out
      location_scores.to_a.sort{|a,b| a[0].downcase<=>b[0].downcase}.each do |array|
        output_line.push array[1]
      end
    end
    wolfer.call('animal',Bio::PSORT::WoLF_PSORT::ANIMAL_LOCATIONS)
    wolfer.call('fungi',Bio::PSORT::WoLF_PSORT::FUNGI_LOCATIONS)
    wolfer.call('plant',Bio::PSORT::WoLF_PSORT::PLANT_LOCATIONS)
    
    #  Plasmit(2): Localisation
    headers.push 'plasmit' if do_headers
    if falciparum_plasmit_predictions[plasmodb].nil?
      $stderr.puts "Warning: No PlasMit prediction found for `#{plasmodb} '"
    end
    output_line.push falciparum_plasmit_predictions[plasmodb]
    
    #  Number of C. hominis Genes in Official Orthomcl Group(1): Localisation
    if do_headers
      headers.push 'Chominis_orthologues'
      headers.push 'Cparvum_orthologues'
      headers.push 'Cmuris_orthologues'
    end
    output_line.push !falciparum_chom_orthologues[plasmodb].nil?
    output_line.push !falciparum_cpar_orthologues[plasmodb].nil?
    output_line.push !falciparum_cmur_orthologues[plasmodb].nil?
    
    #  Chromosome(14): Localisation
    # Encode as 14 dummy variables
    if do_headers
      (1..14).each do |chr| headers.push "chromosome#{chr}"; end
    end
    (1..14).each do |chr|
      output_line.push falciparum_chromosomes[plasmodb].to_s == chr.to_s
    end
    
    
    #  DeRisi 2006 3D7 Timepoint 22(2): Localisation
    #  DeRisi 2006 3D7 Timepoint 23(2): Localisation
    #  DeRisi 2006 3D7 Timepoint 47(2): Localisation
    #  DeRisi 2006 3D7 Timepoint 49(2): Localisation
    interests = [:timepoint22,
       :timepoint23,
       :timepoint47,
       :timepoint49,
       :amplitude]
    if do_headers
      interests.each do |time|
        headers.push time.to_s
      end
      headers.push 'distanceTo48hrs'
    end
    if falciparum_lifecycle_data[plasmodb]
      lifecycle_sample = falciparum_lifecycle_data[plasmodb].sample
      [22,23,47,49].each do |i|
        output_line.push lifecycle_sample[:timepoints][i-1]
      end
      output_line.push lifecycle_sample[:amplitude].to_f
      # Distance to typical invasion maximum hr.
      supplemented_array = [lifecycle_sample[:timepoints],lifecycle_sample[:timepoints]].flatten
      max_hr = nil
      max_average = 0
      (3..55).each do |possible|
        average = supplemented_array[possible-2..possible+2].sum/5.0
        if max_average < average
          hr = possible % 53 #54=> 1 and 55=>2
          max_hr = hr
          max_average = average
        end
      end
      # Now we have the maximum hr, how far away from 48 is it?
      distance_from_max_invasion = nil
      if max_hr < 22 #48-53/2 = 21.5
        distance_from_max_invasion = 53-48+max_hr
      else
        distance_from_max_invasion = (48-max_hr).abs
      end
      output_line.push distance_from_max_invasion
    else
      4.times do output_line.push 1; end #average timepoint
      output_line.push 2.3731476877 #average ampltude
      output_line.push 12 #approx half-way between 48 and max distance to 48
    end
    
    
    # proteomes
    headers.push 'nuclear_proteome' if do_headers
    output_line.push nuclear_proteome_plasmodb_ids.include?(plasmodb)
    
    # number of acidic / basic residue in the first 25 amino acids
    acidics = 0
    protein_sequence[0..24].each_char do |aa|
      acidics += 1 if %w(D E).include?(aa)
    end
    headers.push 'acidics_in_first25' if do_headers
    output_line.push acidics
    
    basics = 0
    protein_sequence[0..24].each_char do |aa|
      basics += 1 if %w(R K H).include?(aa)
    end
    headers.push 'basics_in_first25' if do_headers
    output_line.push basics
    
    # number of transmembrane domains, after the Signal peptide has been removed
    tmhmm = Bio::TMHMM::TmHmmWrapper.new
    result = tmhmm.calculate(signalp.cleave(protein_sequence))
    headers.push 'TMDs' if do_headers
    output_line.push result.transmembrane_domains.length > 0

    # Contained in HP1 associated genomic regions?
    headers.push 'hp1' if do_headers
    output_line.push hp1_positive_genes.include?(plasmodb)
    
    # In un-aligned regions of the genome (sub-telomeric)?
    headers.push 'unconserved_genomic_region' if do_headers
    output_line.push unconserved_region_genes.include?(plasmodb)

    # dcnls
    headers.push 'dcNLS_4of5' if do_headers
    if dcnls4of5[plasmodb] and [2,3,4].include?(dcnls4of5[plasmodb].length)
      output_line.push true
    else
      output_line.push false
    end
    headers.push 'dcNLS_5of6' if do_headers
    if dcnls5of6[plasmodb] and [2].include?(!dcnls5of6[plasmodb].length)
      output_line.push true
    else
      output_line.push false
    end

    # conserved 5' end when blasted against toxo orthologue?
    headers.push 'falciparum_toxo_blast_start' if do_headers
    if falciparum_to_toxo_min_starts[plasmodb]
      output_line.push falciparum_to_toxo_min_starts[plasmodb]
    else
      output_line.push 500 
    end
    
    # conserved between 20 and 100 amino acids after the start?
    headers.push 'falciparum_toxo_blast_start_20_to_100' if do_headers
    if falciparum_to_toxo_min_starts[plasmodb]
      start = falciparum_to_toxo_min_starts[plasmodb]
      output_line.push start>20 and start<100
    else
      output_line.push false
    end
    
    # Output headers
    puts headers.join(',') if do_headers
    do_headers = false
    
    # Convert
    output_line.flatten!
    output2 = []
    output_line.each_with_index do |output, i|
      if output.kind_of?(TrueClass)
        output2.push 1
      elsif output.kind_of?(FalseClass)
        output2.push 0
      elsif output.kind_of?(Numeric)
        output2.push output
      elsif i==0 # Special case for the first column which is the PlasmoDB ID
        output2.push output
      else
        pp output_line
        raise Exception, "Unexpected output class for #{output.inspect} when investigating #{plasmodb}, index #{i}\nTotal output"
      end
    end
    puts output2.join(",")
    
    # progress
    progress.inc
  end
  progress.finish
end
