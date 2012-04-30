#!/usr/bin/env ruby

require 'rubygems'
require 'bio'
require 'optparse'
require 'csv'

options = {
  :email => nil,
}
o = OptionParser.new do |opts|
  opts.banner = "Usage: entrez_proteins_from_nucleotide_ids.rb --email/-e <email_for_ncbi> <nucleotide_accessions_file>\n\nScript to extract a FASTA file of protein sequences associated with input GenBank nucleotide accession numbers. The input file can either be a one column input, or 2 column tab separated with the second column specifying a lable to prepend each translated sequence with."

  opts.on('-e', "--email EMAIL", "Email address to provide to NCBI") do |v|
    options[:email] = v
  end

  opts.on('-g', "--genbank GENBANK_FILE", "Skip entrez, instead provide a GenBank file to extract the CDS from directly. No input file of accessions is used either.") do |v|
    options[:genbank] = v
  end
end.parse!

unless options[:email]
  $stderr.puts "ERROR: email address required\n\n"
  $stderr.puts o.banner
  exit 1
end


process_genbank = lambda do |gb_string, prepending|
  gb = Bio::GenBank.new(gb_string)
  gb.each_cds do |cds|
    ha = cds.to_hash

    # skip pseudogenes
    next if ha['protein_id'].nil?
    
    # define the name as >protein_id description
    name = ""
    unless prepending.nil?
      name = "#{prepending}|"
    end
    name = "#{name}#{ha['protein_id'][0]} #{ha['product'][0]}"

    # Get the protein translation
    translates = ha['translation']
    if translates.length == 1
      puts ">#{name}"
      puts translates[0]
    else
      $stderr.puts "INFO Skipping #{name}, since #{translates.length} 'translation's were found, expected 1"
    end
  end
end



if options[:genbank]
  process_genbank.call(File.open(options[:genbank]).read, 'abc')

else
  CSV.foreach(ARGV[0],:col_sep => "\t") do |row|
    line = row[0].strip
    next if line.length == 0
    prepend = row[1]
    
    $stderr.puts "Fetching #{prepend}|#{line} from GenBank.."
    gb = Bio::NCBI::REST.efetch(line, {"db"=>"nuccore","rettype"=>"gb", "email"=>options[:email]})

    process_genbank.call(gb, prepend)
  end
end
