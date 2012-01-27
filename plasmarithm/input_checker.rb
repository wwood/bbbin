#!/usr/bin/env ruby

# Data-specific tests for Plasmarithm input generation

require 'csv'

if __FILE__ == $0
  tests = {
  'PFL1960w' => { # a var gene
    'signalp' => 0,
    'plasmit' => 0,
    'chom_orthologues' => 0,
    'MPMP_PfEMP1' => 1,
    'MPMP_Plasmodium falciparum histone acetylation and methylation' => 0,
    'timepoint22' => 0.0,
    'distanceTo48hrs' => 12,
    'hp1' => 1,
    'unconserved_genomic_region' => 1,
    'exportpredKLD' => 1,
    'falciparum_toxo_blast_start' => 500,
    
  },
  
  'PFI1475w' => {#MSP1
    'signalp' => 1,
    'aliphatic_index' => 85.73837209302324,
    'gravy' => '-0.6550581395348859',
    'isoelectric_point' => 6.51,
    'distanceTo48hrs' => 0,
    'falciparum_toxo_blast_start_20_to_100' => 0,
    'maurers_cleft' => 1,
    'acidics_in_first25' => 1,
    'basics_in_first25' => 2,
    'acidics_in_first25_afterSPcleavage' => 5,
    'basics_in_first25_afterSPcleavage' => 3,
    'sir2_wild_type_schizont' => 96.5,
    'InterPro_IPR007740' => 0, #Mitochondrial large subunit ribosomal protein (Img2)

      
    },
  'MAL13P1.200' => { #mito ribosomal protein
    'plasmit' => 1,
    'chom_orthologues' => 0,
    'distanceTo48hrs' => 26,
    'hp1' => 0,
    'unconserved_genomic_region' => 0,
    'falciparum_toxo_blast_start_20_to_100' => 1,
    'maurers_cleft' => 0,
    'exportpredKLD' => 0,
    'InterPro_IPR007740' => 1, #Mitochondrial large subunit ribosomal protein (Img2)
    

    }
  }
  
  
  # read output csv
  counted = 0
  CSV.foreach(ARGV[0], :headers => true) do |row|
    if tests[row['plasmodb']]
      tests[row['plasmodb']].each do |tester, answer|
        unless (row[tester].to_s == answer.to_s)
          $stderr.puts "Failure! for #{row['plasmodb']}: #{tester}, expected #{answer}, got #{row[tester]}"
        end
      end
      counted += 1
    end
  end
  $stderr.puts "Unexpected number of genes to test, expected #{tests.length}, found #{counted}" unless counted == tests.length
end