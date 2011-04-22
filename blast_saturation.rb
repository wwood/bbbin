#!/usr/bin/env ruby

# This script is similar to the original NCBI blast. The difference is that it tries to remove 
# problems with one area of the protein hogging all the hits.
#

ACCEPTABLE_BLAST_PLUS_PROGRAMS = %w(blastp blastn tblastn tblastx blastx)

class Blaster
  # Using (if true) BLAST+ or (if false) simply the older blast2 suite of programs?
  attr_reader :blast_plus_program
  attr_reader :verbose
  
  def initialize(blast_plus_program=nil)
    @blast_plus_program = blast_plus_program
    @verbose = true
  end
  
  def hush
    @verbose = false
  end
  
  def blast(fasta_filename, blast_options, output_filename, blast_iteration_number)
    $stderr.puts "Running blast ##{blast_iteration_number}, starting at #{`date`}" if @verbose
    # Do the actual blast
    cmd = nil
    if @blast_plus_program.nil?
      cmd = "blastall -m 8 -i '#{fasta_filename}' -o #{output_filename} #{blast_options}"
    else
      cmd = "#{@blast_plus_program} -outfmt 6 -query '#{fasta_filename}' -out #{output_filename} #{blast_options}"
    end
    $stderr.puts "BLAST CMD LINE$ #{cmd}" if @verbose
    `#{cmd}`
    print `cat #{output_filename}` #output to the user each of the blast results
    $stderr.puts "Blast ##{blast_iteration_number} returned #{`wc -l #{output_filename}`.to_i} result lines" if @verbose
  end
end



if __FILE__ == $0
  require 'tempfile'
  
  # Parse options
  USAGE = [
  'Usage: blast_saturation.rb [-q] [-p <blast+_program>] -f <fasta_filename> -b <blast_options>',
  
  ]
  options = {
    :blast_plus_program => nil,
    :verbose => true
  }
  require 'optparse'
  o = OptionParser.new do |opts|
    opts.banner = USAGE
    
    opts.on('-f','--fasta-filename FASTA_FILENAME', "fasta filename of sequences to BLAST") do |v|
      options[:fasta_filename] = v
    end    
    
    verbose = "Options to pass to the blast program. For BLAST+ usage (-p), required argument is -db, illegal arguments that give undefined behaviour are -out, -outfmt or -query. If using the older version of blast, -p and -d are required, while -m, -o and -i are illegal."
    opts.on('-b','--blast-options BLAST_OPTIONS', verbose) do |v|
      #TODO: make the illegal things really illegal by catching them
      options[:blast_options] = v
    end
    
    opts.on("-p", "--blast_plus_program [PROGRAM]", "Use the specified blast+ program instead of the previous BLAST2 version.") do |v|
      unless ACCEPTABLE_BLAST_PLUS_PROGRAMS.include?(v)
        $stderr.puts "Unexpected blast+ program found: `#{v}'. Expected one of #{ACCEPTABLE_BLAST_PLUS_PROGRAMS.sort.join(", ")}."
        exit 1
      end
      options[:blast_plus_program] = v
    end
    
    opts.on('-q','--quiet','Opposite of verbose. Default is not quiet (verbose is on)') do
      options[:verbose] = false
    end
  end
  o.parse!
  
  # require fasta input
  # should be able to set arbitrary blast parameters (required, otherwise blast will not work - db and program at least needed).
  unless options[:fasta_filename] and options[:blast_options]
    $stderr.puts USAGE
    exit 1
  end
  fasta_filename = options[:fasta_filename]
  blast_options = options[:blast_options]
  
  # setup blaster, which is differently done when blast+ has been specified
  blaster = Blaster.new options[:blast_plus_program]
  blaster.hush unless options[:verbose]
  
  # should be able to set number deep at each point #TODO
  
  blast_iteration_number = 1
  # run the initial blast
  blast_result_file = Tempfile.new("blast_saturation#{blast_iteration_number}")
  blast_result_file.close
  
  blaster.blast(fasta_filename, blast_options, blast_result_file.path, blast_iteration_number)
  
  last_masked_file = File.new(fasta_filename,'r')
  # while there are proteins still being hit
  while File.size(blast_result_file.path) > 0
    # mask out the hit areas
    $stderr.puts "Masking the fasta file parts that have blast hits, starting at #{`date`}" if options[:verbose]
    masked_fasta_file = Tempfile.new("masked_fasta#{blast_iteration_number}")
    `blast_mask.rb -m #{last_masked_file.path} #{blast_result_file.path} >#{masked_fasta_file.path}`
    #$stderr.puts `cat #{masked_fasta_file.path}` #debug
    blast_iteration_number += 1
    # redo the blast, this time using the masked sequence as input
    blast_result_file = Tempfile.new("blast_saturation#{blast_iteration_number}")
    
    # do the actual blasting
    blaster.blast(masked_fasta_file.path, blast_options, blast_result_file.path, blast_iteration_number)
    
    # Get ready for next time. It's the old 'last = cur' line - just like programming linked lists in ugrad
    last_masked_file = masked_fasta_file
  end
end