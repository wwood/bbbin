#!/usr/bin/env ruby

# This script is similar to blast. The difference is that it tries to remove 
# problems with one area of the protein hogging all the hits.
#


class Blaster
  def blast(fasta_filename, blast_options, output_filename, blast_iteration_number)
    $stderr.puts "Running blast ##{blast_iteration_number}, starting at #{`date`}"
    # Do the actual blast
    cmd = "blastall -m 8 -i '#{fasta_filename}' -o #{output_filename} #{blast_options}"
    `#{cmd}`
    print `cat #{output_filename}` #output to the user each of the blast results
    $stderr.puts "Blast ##{blast_iteration_number} returned #{`wc -l #{output_filename}`.to_i} result lines"
  end
end



if __FILE__ == $0
  require 'tempfile'
  
  # Parse options
  USAGE = [
  'Usage: blast_saturation.rb <fasta_filename> <blast_options>',
  '  where <fasta_filename> is the fasta file of sequences to blast, and <blast_options> is a string you might give to blastall.',
  '   <blast_options> must include -d and -p, but cannot include -m or -i'
  ]
  #  options = {
  #  :print_masked_sequences_only => false,
  #  }
  #  o = OptionParser.new do |opts|
  #    opts.banner = USAGE
  #    
  #    opts.on("-m", "--masked-only", "Print out only those sequences that have blast hits (or equivalently, only those that are masked") do |v|
  #      options[:print_masked_sequences_only] = true
  #    end
  #  end
  #  o.parse!
  
  # require fasta input
  # should be able to set arbitrary blast parameters (required, otherwise blast will not work - db and program at least needed).
  unless ARGV.length == 2
    $stderr.puts USAGE
    exit 1
  end
  fasta_filename = ARGV[0]
  blast_options = ARGV[1]
  blaster = Blaster.new
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
    $stderr.puts "Masking the fasta file parts that have blast hits, starting at #{`date`}"
    masked_fasta_file = Tempfile.new("masked_fasta#{blast_iteration_number}")
    `blast_mask.rb -m #{last_masked_file.path} #{blast_result_file.path} >#{masked_fasta_file.path}`
    blast_iteration_number += 1
    # redo the blast, this time using the masked sequence as input
    blast_result_file = Tempfile.new("blast_saturation#{blast_iteration_number}")
    
    # do the actual blasting
    blaster.blast(masked_fasta_file.path, blast_options, blast_result_file.path, blast_iteration_number)
    
    # Get ready for next time. It's the old 'last = cur' line - just like programming linked lists in ugrad
    last_masked_file = masked_fasta_file
  end
end