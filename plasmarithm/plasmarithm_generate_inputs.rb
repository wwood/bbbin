#!/usr/bin/env ruby

# biology related gems
require 'bio'
require 'bio-plasmoap'
require 'bio-signalp'
require 'reubypathdb'
require 'bio-tm_hmm'
require 'bio-exportpred'
require 'bio-wolf_psort_wrapper'
require 'bio-isoelectric_point'
require 'bio-aliphatic_index'
require 'bio-hydropathy'
# regular gems
require 'progressbar'
# stdlib includes
require 'pp'


# Generated using plasmit_screenscraper.rb
PLASMIT_PREDICTIONS_FILE = '/home/ben/phd/data/Plasmodium falciparum/plasmit/PlasmoDB8.2.plasmit_screenscraped.csv'
# Generated using plasmit_screenscraper.rb, after removing the signal peptides
PLASMIT_SP_CLEAVED_PREDICTIONS_FILE = '/home/ben/phd/data/Plasmodium falciparum/plasmit/PlasmoDB8.2.signalp_cleaved.plasmit_screenscraped.csv'
# Generated using orthomcl_species_jumper.rb
ORTHOLOGUES_FILE = '/home/ben/phd/data/Plasmodium falciparum/genome/PlasmoDB/8.2/orthologues.txt'
# Downloaded directly PlasmoDB
FALCIPARUM_GENE_INFORMATION_PATH = '/home/ben/phd/data/Plasmodium falciparum/genome/PlasmoDB/8.2/PfalciparumGene_PlasmoDB-8.2.txt'
# Supplementary data to the 2006 microarray paper
DERISI_3D7_MICROARRAY_LIFECYCLE_PATH = '/home/ben/phd/data/Plasmodium falciparum/microarray/DeRisi2006/3D7_overview_S6.csv'
# Obtained from Till Voss
HP1_POSITIVE_GENES_PATH = '/home/ben/phd/data/Plasmodium falciparum/HP1chip/hp1.txt'
# Obtained by running falciparum_unconserved_regions.rb
FALCIPARUM_UNCONSERVED_REGION_PLASMODB_IDS = 'falciparum_unconserved_region_genes.window5max1.csv'
# used dcnls.rb and some other things for these two
DCNLS_4_OF_5_PATH = 'dcnls.4of5.csv'
DCNLS_5_OF_6_PATH = 'dcnls.5of6.csv'
# Conservation of 5' end of the gene
FALCIPARUM_TO_TOXO_BLAST_CSV = 'falciparum8.2versusToxo8.2.blast.csv'
# Sir2 A and B microarray data - Downloaded from PlasmoDB 8.2 with New Search>Search for Genes>Transcript Expression>Microarray evidence
# then selecting "P.f. sir2 KO (percentile)", min expression 1, max expression 100, and then downloading the table as a CSV, not downloading the "Product Description" column
# For WT ones, chose min 0, max 100 instead.
SIR2_KO_MICROARRAYS = {
:sir2a_ring => 'input_data/Sir2AringPercentiles.csv',
:sir2a_trophozoite => 'input_data/Sir2AtrophozoitePercentiles.csv',
:sir2a_schizont => 'input_data/Sir2AschizontPercentiles.csv',
:sir2b_ring => 'input_data/Sir2BringPercentiles.csv',
:sir2b_trophozoite => 'input_data/Sir2BtrophozoitePercentiles.csv',
:sir2b_schizont => 'input_data/Sir2BschizontPercentiles.csv',
:sir2_wild_type_ring => 'input_data/Sir2wildTypeRingPercentiles.csv',
:sir2_wild_type_trophozoite => 'input_data/Sir2wildTypeTrophozoitePercentiles.csv',
:sir2_wild_type_schizont => 'input_data/Sir2wildTypeSchizontPercentiles.csv',
}
# Invasion pathway KO data, derived from PlasmoDB 8.2 similarly to above, but instead using a 2 fold up/down regulation cutoff, compared to the reference experiment "3D7 WT 48hr"
INVASION_KO_MICROARRAYS = {
:eba140down => 'input_data/eba140downregulated2fold.txt',
:eba140up => 'input_data/eba140upregulated2fold.txt',
:eba175down => 'input_data/eba175downregulated2fold.txt',
:eba175up => 'input_data/eba175upregulated2fold.txt',
:rh2b_down => 'input_data/rh2bDownregulated2fold.txt',
:rh2b_up => 'input_data/rh2bUpregulated2fold.txt',
}
PROTEOMES = {
  # Oehring and Woodcroft et. al. 2012 (hopefully - not yet accepted, but on PlasmoDB at least)
:core_nuclear => '/home/ben/phd/voss_proteome/June2010/nuclear_bioinfo_final6.ids',
  # Modified by taking the PlasmoDBs from the manuscript itself.
:food_vacuole => 'input_data/food_vacuole_proteome.txt',
  # Modified from table 2
:maurers_cleft => '/home/ben/phd/data/Plasmodium falciparum/proteomics/MaurersCleft2005/table2.ids',
}
# Blast hits. Generated using blast_inputs_creation.rb
pros = %w(
PbergheiAnnotatedProteins_PlasmoDB-8.2.fasta
PvivaxAnnotatedProteins_PlasmoDB-8.2.fasta

TgondiiME49AnnotatedProteins_ToxoDB-7.2.fasta

CmurisAnnotatedProteins_CryptoDB-4.6.fasta
ChominisAnnotatedProteins_CryptoDB-4.6.fasta
CparvumAnnotatedProteins_CryptoDB-4.6.fasta

BbovisT2BoAnnotatedProteins_PiroplasmaDB-1.1.fasta
TannulataAnkaraAnnotatedProteins_PiroplasmaDB-1.1.fasta
TparvaMugugaAnnotatedProteins_PiroplasmaDB-1.1.fasta
)
BLAST_RESULT_FILENAMES = {}
%w(pberghei pvivax tgondii cmuris chominis cparvum bbovis tannulata tparva).each_with_index do |sp, i|
  pro = pros[i]
  BLAST_RESULT_FILENAMES[sp] = "input_data/falciparum8.2versus#{pro}.blast.csv" 
end


class Array; def sum; self.inject{|sum,x| sum + x }; end; end



if __FILE__ == $0
  # Take a file that contains 1 PlasmoDB ID per line, and generate a matrix that 
  # corresponds to the inputs for use in the R package glmnet.
  
  USAGE = 'plasmarithm_generate_inputs.rb <plasmodb_id_list_file>'
  
  
  # Cache blast hits
  blast_hits = {} # Hash of target species => plasmodbid => true
  # e.g. psu|MAL13P1.200 Genbank|BBOV_III001790  41.41 99  57  1 14  111 5 103 3e-18 84.3
  all_blast_target_species = BLAST_RESULT_FILENAMES.keys.sort
  all_blast_target_species.each do |sp|
    File.open(BLAST_RESULT_FILENAMES[sp]).each_line do |line|
      first = line.strip.split("\t")[0]
      matches = first.match(/^psu\|([^\t]+)/)
      raise Exception, "Unexpected blast line found: #{line} / #{first}" unless matches
      plasmodb = matches[1]
      blast_hits[sp] ||= {}
      blast_hits[sp][plasmodb] = true
    end
  end
  all_blast_target_species.each do |sp|
    $stderr.puts "Cached #{blast_hits[sp].length} straight out BLAST hits to #{sp}"
  end
  

  
  # Cache Sir2 KO percentiles
  sir2_percentiles = {} #hash of experiment name => plasmodb => percentile
  SIR2_KO_MICROARRAYS.each do |experiment_name, filename|
    sir2_percentiles[experiment_name.to_s] = {}
    File.open(filename).readlines.each do |l|
      splits = l.split("\t")
      plasmodb = splits[0]
      next if plasmodb == '[Gene ID]'
      percentile = splits[1].to_f
      sir2_percentiles[experiment_name.to_s][plasmodb] = percentile
    end
  end
  all_sir2_microarrays = sir2_percentiles.keys.sort
  
  # Cache Invasion KO 2fold enrichment lists
  invasion_microarrays = {}
  INVASION_KO_MICROARRAYS.each do |experiment_name, filename|
    invasion_microarrays[experiment_name.to_s] = {}
    File.open(filename).readlines.each do |l|
      plasmodb = l.strip
      next if plasmodb == '[Gene ID]'
      invasion_microarrays[experiment_name.to_s][plasmodb] = true
    end
  end
  all_invasion_microarrays = invasion_microarrays.keys.sort
  
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
      raise Exception, "Unexpected PlasMit output: `#{splits[1]}'"
    end
    falciparum_plasmit_predictions[gene_id] = classification
  end
  $stderr.puts "Cached #{falciparum_plasmit_predictions.length} PlasMit predictions"
  
  # Cache Crypto orthologues file
  falciparum_orthologues = {} # hash of target non-falciparum species => plasmodb => output chunk
  file_counter = 0
  all_orthologue_species = [:chom, :cpar, :cmur, :tgon, :bbov, :tann, :tpar, :pber, :pviv]
  File.open(ORTHOLOGUES_FILE).each_line do |line|
    splits = line.strip.split("\t")
    raise unless splits.length >= 2 or splits.length > 5
    plasmodb = splits[2]
    raise if falciparum_orthologues[plasmodb]
    all_orthologue_species.each_with_index do |sp, i|
      index = 2+i
      falciparum_orthologues[sp] ||= {}
      falciparum_orthologues[sp][plasmodb] = splits[index]
    end
    file_counter += 1
  end
  $stderr.puts "Cached sets of orthologues for #{file_counter} P. falciparum sequences"
  
  # Cache chromosome numbers for each gene
  falciparum_chromosomes = {}
  metabolic_pathways = {} # hash of PlasmoDB ID => array of metabolic pathways
  gene_starts = {}
  number_of_exons = {}
  interpro_domains = {}
  EuPathDBGeneInformationTable.new(File.open(FALCIPARUM_GENE_INFORMATION_PATH)).each do |info|
    plasmodb = info.get_info('Gene ID')
    chromosome = info.get_info('Chromosome')
    raise Exception, "found PlasmoDB ID `#{plasmodb}' multiple times in the Gene information file" if falciparum_chromosomes[plasmodb]
    falciparum_chromosomes[plasmodb] = chromosome.to_i
    
    # MPMP
    metabolic_pathways[plasmodb] = info.get_table('Metabolic Pathways').collect{|t| "MPMP_#{t['pathway_id']}"}
    
    # Gene position (for in relation to the chromosome ends (min start of exon)
    gene_starts[plasmodb] = info.get_table('Gene Model').select{|entry|
      entry['Type']=='exon'
    }.collect{|entry|
      entry['Start']
    }.min.to_i
    
    # Number of exons
    number_of_exons[plasmodb] = info.get_table('Gene Model').select{|entry|
      entry['Type']=='exon'
    }.length
    
    # Interpro domains
    interpro_domains[plasmodb] = info.get_table('InterPro Domains').collect{|entry|
      entry['Interpro ID']
    }.uniq
  end
  all_metabolic_pathways = metabolic_pathways.values.flatten.uniq.sort
  all_interpro_domains = interpro_domains.values.flatten.uniq.sort.reject{|i| i.to_s == ''}
  $stderr.puts "Cached #{falciparum_chromosomes.length} chromosome numbers for P. falciparum genes"
  $stderr.puts "Cached #{metabolic_pathways.length} proteins with MPMP"
  
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
  
  # Cache proteomes
  proteome_entrants = {} # hash of proteome name to array of plasmodb IDs
  PROTEOMES.each do |name, filename|
    raise if proteome_entrants[name]
    proteome_entrants[name] = []
    File.open(filename).each_line do |line|
      pl = line.strip
      proteome_entrants[name].push pl unless pl == ''
    end
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
  
  ##### Cache the length of the chromosomes
  chromosome_lengths = {} # hash of chromosome name to length
  # Copied this from the top of the GFF file manually.
  top_of_gff = <<END_OF_TOP
##sequence-region psu|Pf3D7_01  1 643292
##sequence-region psu|Pf3D7_02  1 947102
##sequence-region psu|Pf3D7_03  1 1060087
##sequence-region psu|Pf3D7_04  1 1204112
##sequence-region psu|Pf3D7_05  1 1343552
##sequence-region psu|Pf3D7_06  1 1418244
##sequence-region psu|Pf3D7_07  1 1501717
##sequence-region psu|Pf3D7_08  1 1419563
##sequence-region psu|Pf3D7_09  1 1541723
##sequence-region psu|Pf3D7_10  1 1687655
##sequence-region psu|Pf3D7_11  1 2038337
##sequence-region psu|Pf3D7_12  1 2271478
##sequence-region psu|Pf3D7_13  1 2895605
##sequence-region psu|Pf3D7_14  1 3291871
END_OF_TOP
  top_of_gff.split("\n").each do |line|
    splits = line.strip.split ' '
    regex = /^psu\|Pf3D7_(\d\d)$/
    raise unless splits[1].match(regex)
    chromosome_number = splits[1].match(regex)[1].to_i
    chromosome_lengths[chromosome_number] = splits[3].to_i
    raise Exception, "Found unexpectedly low chromosome length for `#{chromosome_number}' from #{splits[1]}: #{chromosome_lengths[chromosome_number]}" if chromosome_lengths[chromosome_number] < 400
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
    signalp_cleaved = signalp.cleave(protein_sequence)
    
    #  PlasmoAP Score(2): Localisation
    # We already know whether there is a SignalP prediction or not, so just go with that. Re-calculating takes time.
    headers.push 'plasmoap' if do_headers
    output_line.push Bio::PlasmoAP.new.calculate_score(protein_sequence, signalp.signal?, signalp_cleaved).points
    
    #  ExportPred?(2): Localisation
    headers.push ['exportpredRLE', 'exportpredKLD'] if do_headers
    output_line.push Bio::ExportPred::Wrapper.new.calculate(protein_sequence, :no_KLD => true).predicted_rle?
    output_line.push Bio::ExportPred::Wrapper.new.calculate(protein_sequence, :no_RLE => true).predicted_kld?
    
    headers.push 'aliphatic_index' if do_headers
    output_line.push Bio::Sequence::AA.new(protein_sequence).aliphatic_index
    
    headers.push 'gravy' if do_headers
    output_line.push Bio::Sequence::AA.new(protein_sequence).gravy
    
    #  WoLF_PSORT 
    # Encode each with the scores that they were given, accounting for dual locations
    wolfer = lambda do |lineage, all_locations, header_pre, sequence|
      all_locations.collect{|l| l.downcase}.sort.each{|l| headers.push "#{header_pre}wolf_#{lineage}_#{l}"} if do_headers
      result = Bio::PSORT::WoLF_PSORT::Wrapper.new.run(sequence, lineage)
      
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
#    wolfer.call('animal',Bio::PSORT::WoLF_PSORT::ANIMAL_LOCATIONS, '', protein_sequence)
#    wolfer.call('fungi',Bio::PSORT::WoLF_PSORT::FUNGI_LOCATIONS, '', protein_sequence)
#    wolfer.call('plant',Bio::PSORT::WoLF_PSORT::PLANT_LOCATIONS, '', protein_sequence)
#    wolfer.call('animal',Bio::PSORT::WoLF_PSORT::ANIMAL_LOCATIONS, 'SPcleaved_', signalp_cleaved)
#    wolfer.call('fungi',Bio::PSORT::WoLF_PSORT::FUNGI_LOCATIONS, 'SPcleaved_', signalp_cleaved)
#    wolfer.call('plant',Bio::PSORT::WoLF_PSORT::PLANT_LOCATIONS, 'SPcleaved_', signalp_cleaved)
    
    #  PlasMit
    headers.push 'plasmit' if do_headers
    if falciparum_plasmit_predictions[plasmodb].nil?
      $stderr.puts "Warning: No PlasMit prediction found for `#{plasmodb} '"
    end
    output_line.push falciparum_plasmit_predictions[plasmodb]
    
    #  Number Genes in Official OrthoMC Group for various species
    all_orthologue_species.each do |sp|      
      if do_headers
        headers.push "#{sp}_orthologues"
      end
      output_line.push !falciparum_orthologues[sp][plasmodb].nil?
    end
    
    #  Chromosome(14): Localisation
    ## Encode as 14 dummy variables. Not doing this any more because we don't believe they are predictive and are just noise.
    #    if do_headers
    #     (1..14).each do |chr| headers.push "chromosome#{chr}"; end
    #    end
    #     (1..14).each do |chr|
    #      output_line.push falciparum_chromosomes[plasmodb].to_s == chr.to_s
    #    end
    
    # Distance from chromosome ends
    min_distance = gene_starts[plasmodb]
    raise "No chromosome found for #{plasmodb}" if falciparum_chromosomes[plasmodb].nil?
    distance_from_arbitrary_end = chromosome_lengths[falciparum_chromosomes[plasmodb]]-gene_starts[plasmodb]
    min_distance = distance_from_arbitrary_end if distance_from_arbitrary_end < min_distance
    [100000,200000].each do |cutoff|
      kb = cutoff/1000
      headers.push "within#{kb}_kb_of_chromosome_end" if do_headers
      output_line.push min_distance < cutoff
    end
    
    # number of exons
    headers.push 'num_exons' if do_headers
    exons = number_of_exons[plasmodb]
    if exons.nil?
      output_line.push 1
    else
      output_line.push exons
    end
    headers.push '2exons' if do_headers
    output_line.push number_of_exons[plasmodb]==2
    
    # isoelectric point
    headers.push 'isoelectric_point' if do_headers
    output_line.push Bio::Sequence::AA.new(protein_sequence).isoelectric_point
    
    #  Bozdech/DeRisi 2006 Microarray data
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
        output_line.push Math.log(lifecycle_sample[:timepoints][i-1])
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
      4.times do output_line.push Math.log(1); end #average timepoint
      output_line.push 2.3731476877 #average ampltude
      output_line.push 12 #approx half-way between 48 and max distance to 48
    end
    
    
    # proteomes
    proteome_entrants.sort{|a,b| a[0]<=> b[0]}.each do |arr|
      name = arr[0]
      plasmodb_ids = arr[1]
      headers.push name if do_headers
      output_line.push plasmodb_ids.include?(plasmodb)
    end
    
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
    
    # Number of acidics/basics after the SP has been removed
    acidics = 0
    signalp_cleaved[0..24].each_char do |aa|
      acidics += 1 if %w(D E).include?(aa)
    end
    headers.push 'acidics_in_first25_afterSPcleavage' if do_headers
    output_line.push acidics
    
    basics = 0
    signalp_cleaved[0..24].each_char do |aa|
      basics += 1 if %w(R K H).include?(aa)
    end
    headers.push 'basics_in_first25_afterSPcleavage' if do_headers
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
      output_line.push (start>20 and start<100)
    else
      output_line.push false
    end
    
    # Sir2 arrays
    all_sir2_microarrays.each do |experiment_name|
      headers.push experiment_name if do_headers
      if sir2_percentiles[experiment_name][plasmodb]
        output_line.push sir2_percentiles[experiment_name][plasmodb]
      else
        output_line.push 50 #hopefully not very prevalent
      end
    end
    
    # Invasion arrays
    all_invasion_microarrays.each do |experiment_name|
      headers.push experiment_name if do_headers
      output_line.push !(invasion_microarrays[experiment_name][plasmodb].nil?)
    end
    
    # blast hits
    all_blast_target_species.each do |sp|
      headers.push "blast_#{sp}" if do_headers
      output_line.push blast_hits[sp].key?(plasmodb)
    end
    
    # Metabolic pathways
    all_metabolic_pathways.each do |mpmp| 
      headers.push mpmp if do_headers
      output_line.push (metabolic_pathways[plasmodb] and metabolic_pathways[plasmodb].include?(mpmp)) 
    end
    
    # InterPro domains
    all_interpro_domains.each do |ipr|
      headers.push "InterPro_#{ipr}" if do_headers
      output_line.push interpro_domains[plasmodb].include?(ipr)
    end
    
    
    
    
    ##################################################################
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
        raise Exception, "Unexpected output class for #{output.inspect} when investigating #{plasmodb}, index #{i} (header? #{headers[i]})\nTotal output"
      end
    end
    puts output2.join(",")
    
    # progress
    progress.inc
  end
  progress.finish
end
