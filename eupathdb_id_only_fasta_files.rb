#!/usr/bin/env ruby


# Output a single file that only has in it the EuPathDB ID, and not all the other stuff that usually comes in the id line of the fasta files

require 'rubygems'
require 'reubypathdb'

DATABASES = %w(PlasmoDB ToxoDB CryptoDB PiroplasmaDB)
    
if __FILE__ == $0
  DATABASES.each do |db, version|
    EuPathDBSpeciesData.species_data_from_database(db).each do |species_data|
      EuPathDBSpeciesData.new(species_data.name, '/home/ben/phd/data').protein_fasta_file_iterator.each do |s|
        print '>'
        puts s.gene_id
        puts s.sequence
      end
    end
  end
end