#!/usr/bin/env ruby

# A script to make the design primers for T. gondii ku80 3' replacement.
# Currently
# this script only works by taking a sequence, and it spits out a list of
# possible enzyme digestion sites. This can be helpful in itself. It then
# tries to find primer sites, but I wouldn't trust those results just yet, and
# it isn't likely to find any suitable primers anyway.
#
# The list of enzymes that are compatible can be changed by modifying this
# script,
# or using the -e option. By default, it is overly ben-specific.

require 'rubygems'
require 'bio'
require 'optparse'
require 'eupathdb_gene_information_table' #from the reubypathdb rubygem

require 'tempfile'
require 'primer_design_extras' # a loose collection of helpful

# classes for primer design

# Pre-requisites
# * a working install of emboss, and in particular the 'restrict' program should
# be runnable and should automatically be able to find the rebase database of
# enzyme restriction sites.
# * common unix-utils
# * an accessible copy of the enzyme_matcher.rb program (from bbbin)
# * the oligotm program, and primer3_core

# Delete this class later, when it is incorporated into a more suitable rubygem
# A class for extracting gene info from a gene info file
class EuPathDBGeneInformationFileExtractor
  def initialize(filename = nil)
    @filename = filename
  end

  # Returns a EuPathDBGeneInformation object corresponding to the wanted key. If
  # there are multiple in the file, only the first is returned.
  def extract_gene_info(wanted_gene_id)
    EuPathDBGeneInformationTable.new(File.open(@filename)).each do |gene|
      return gene if wanted_gene_id == gene.get_info('Gene Id')
    end
    return nil
  end
end

# If we are running this file as a script
if $0 == __FILE__
  options = {
    :gene_information_table_path => '/home/ben/phd/data/Toxoplasma gondii/ToxoDB/6.3/TgondiiME49Gene_ToxoDB-6.3.txt',
    :upstream_sequence => nil,
    :possible_enzymes_to_cut_with_filename => '/home/ben/phd/toxo_experimental_localisations/vectors/LIC HA3:DHFR non-cutters.txt',
    :enzymes_in_freezer_filename => '/home/ben/phd/ralphlab/enzymeList.txt',
    :enzymes_to_avoid => %w(MauBI),
    :print_vector_compatible_sites => false,
    :print_freezer_compatible_sites => false,
    #    :recombination_sequence_buffer => 350, #how far way from priming sites
    # can
    # the restriction site be?
    :five_prime_homologous_recombination_minimum => 800,
    :three_prime_homologous_recombination_minimum => 350,
    :enzyme => nil,
    :primer_temperature_minimum => 50,
    :primer_temperature_optimal => 55,
    :primer_temperature_maximum => 65,
    :primer_gc_clamp => 1, #start out trying to get primers with this gc clamp.
    # If none can be found, decrement until options[:real_gc_clamp_minimum] is
    # reached
    :real_gc_clamp_minimum => 1, #don't try getting primers with anything less
    # than this
    :enzymes_on_order => [], #extra enzymes that will be in the freezer soon,
    :primer_search_area => 500, # Need to choose some finite length to choose
    # primers in, so incidence of non-unique restriction sites is reduced
  }
  wanted_gene_id = nil
  o = OptionParser.new do |opts|
    opts.banner = [
      'Usage: program -s <upstream_sequence>'
    ]
    opts.on("-s", "--sequence FASTA_FILENAME", "A fasta file of the sequence to design primers to (skip the step of extracting these based on the ToxoDB ID)") do |f|
      Bio::FlatFile.open(File.open(f)).each do |e|
        unless options[:upstream_sequence].nil?
          raise Exception, "Only one sequence can be worked on at a time - error in fasta file"
        end
        options[:upstream_sequence] = e.seq.to_s
      end
    end
    opts.on(nil, "--print-vector-compatible", "Print out vector- (but perhaps not freezer-) compatible unique restriction sites") do
      options[:print_vector_compatible_sites] = true
    end
    opts.on(nil, "--print-freezer-compatible", "Print out freezer- (and vector-) compatible unique restriction sites") do
      options[:print_freezer_compatible_sites] = true
    end
    opts.on('-e', "--enzyme ENZYME_NAME", "Use this enzyme to do the cutting") do |e|
      options[:enzyme] = e
    end
    opts.on('-o', "--enzymes_on_order ENZYME_NAMES", "extra enzymes that will be in the freezer soon") do |e|
      options[:enzymes_on_order] = e.split(',')
    end
  end
  o.parse!

  if ARGV.length != 0
    $stderr.puts o.help
    exit
  end

  # ============================================================================
  # unimplemented code to extract the sequence automatically
  upstream_sequence = options[:upstream_sequence]
  if upstream_sequence.nil? #this should neveer
    raise Exception, "not implemented - a fasta file must be specified"

    wanted_gene_id = ARGV[0]

    # Input the position that is to be inserted. Presumably a ToxoDB Gene id,
    # assuming a 3' insertion at the end of the gene just before the stop codon
    # First, extract gene info
    gene_info = EuPathDBGeneInformationFileExtractor.new(options[:gene_information_table_path]).extract_gene_info(wanted_gene_id)

    # Parse +/-ve direction of gene
    location = gene_info.get_info('Genomic Location')
    $stderr.puts "Genomic location: #{location}"
    direction = nil #true => fwd, false => backwards
    chromosome = nil
    if matches = location.match(/(.*): [\d,]+ - [\d,]+ \(([+-])\)$/)
    chromosome = matches[1]
    direction = matches[2]
    else
      raise Exception, "Couldn't parse genomic location line `#{location}'"
    end
    $stderr.puts "Found a #{direction} direction sequence on #{chromosome}"

  # Parse out the end exon - highest ending if +ve direction, lowest ending if -ve direction
  # It doesn't appear to be easily possible to extract the cds positions from the info file, only the transcript exons. Need a separate table?

  # Extract 1kb upstream of insertion point
  end

  $stderr.puts "Found a #{upstream_sequence.length}bp sequence found to match against"

  vector_compatible_cuts = []
  freezer_compatible_cuts = []

  # method to take a file that was output from reas input and return an array
  enzyme_cuts_array_from_filename = lambda do |filename|
    to_return = []
    File.foreach(filename) do |line|
      splits = line.split(' ')
      raise if splits.length != 9
      cut = EnzymeCut.new
      cut.start = splits[5].to_i
      cut.stop = splits[6].to_i
      cut.enzyme = splits[3]
      to_return.push cut
    end
    return to_return
  end

  # ============================================================================
  # Find restriction sites in the sequence, then afterwards parse them into so
  # that:
  # * They are vector-compatible
  # * (maybe) we have the enzymes in the freezer
  Tempfile.open('first_fasta') do |first_fasta_file|
    first_fasta_file.puts '>seq'
    first_fasta_file.puts upstream_sequence
    first_fasta_file.close
    Tempfile.open('restrict') do |restrict_output_file|
    # find the unique restriction sites
      cmd = "restrict -sequence #{first_fasta_file.path} -sitelen 4 -enzymes all -outfile #{restrict_output_file.path}"
      puts `#{cmd}`

      # extract the list of enzyme cutting sites
      Tempfile.open('awk_output') do |awk_out|
        cmd = "awk 'NF==9{print $0}' #{restrict_output_file.path} >#{awk_out.path}"
        `#{cmd}`
        print 'Found this many sites in the whole sequence given, before culling at all: '
        puts `wc -l #{awk_out.path} |awk '{print $1}'`

        # merge the list of enzymes with the ones that are known to not be
        # a) commercially available
        # b) don't also cut the LIC Ku80 vector
        Tempfile.open('enzyme_match1') do |enzyme_match1|
          `enzyme_matcher.rb -f '#{options[:possible_enzymes_to_cut_with_filename]}' #{awk_out.path} >#{enzyme_match1.path}`
          print 'Found this many vector-compatible: '
          puts `wc -l #{enzyme_match1.path} |awk '{print $1}'`
          puts `cat #{enzyme_match1.path}` if options[:print_vector_compatible_sites]

          # parse restriction sites
          vector_compatible_cuts = enzyme_cuts_array_from_filename.call(enzyme_match1.path)

          # Do we have these enzymes on the lab's list of enzymes in the freezer?
          Tempfile.open('enzyme_match2') do |enzyme_match2|
            args = ''
            unless options[:enzymes_on_order].nil? or options[:enzymes_on_order].empty?
              args += "-e #{options[:enzymes_on_order].join(',')}"
            end
            cmd = "enzyme_matcher.rb -f '#{options[:enzymes_in_freezer_filename]}' #{args} #{enzyme_match1.path} >#{enzyme_match2.path}"
            `#{cmd}`
            print 'Found this many freezer-compatible: '
            puts `cat #{enzyme_match2.path} |wc -l`
            puts `cat #{enzyme_match2.path}` if options[:print_freezer_compatible_sites]

            # parse the restriction sites
            freezer_compatible_cuts = enzyme_cuts_array_from_filename.call(enzyme_match2.path)
          end
        end
      end
    end
  end

  # ============================================================================
  # Remove unsuitable enzymes based on their cut position(s)
  #
  # Are unique restriction sites within restricition site-800bp-500bp to 3' end
  # of
  # sequence?
  # If no, then they are of no use, discard that restriction site.
  # This lambda takes an array of EnzymeCut objects, and returns those that are
  # suitable
  remove_unsuitable_enzymes = lambda do |enzyme_cut_array|
    possibles = PossiblyRestrictedNucleotideSequence.new(enzyme_cut_array)
    enzyme_cut_array.select do |c|
      beginner = c.start-options[:five_prime_homologous_recombination_minimum]-options[:primer_search_area]
      if beginner <= 0 #this hopefully won't happen much for restriction sites
      # that
      # end up being useful
      false
      elsif c.start > upstream_sequence.length - options[:three_prime_homologous_recombination_minimum]
      # Enzyme cut sites are too close to the 3' end. Need enough homolgous
      # sequence.
      false
      else
      # remove from array if the enzyme cuts in more than one place.
      answer = possibles.unique_within_region?(c, beginner, upstream_sequence.length)
      answer
      end
    end
  end

  # Anything suitable in the freezer?
  freezer_compatible_cuts = remove_unsuitable_enzymes.call(freezer_compatible_cuts)
  $stderr.puts "After culling for uniqueness, #{freezer_compatible_cuts.length} freezer-compatible sequences remain"
  # Anything in the whole wide world?
  vector_compatible_cuts = remove_unsuitable_enzymes.call(vector_compatible_cuts)
  $stderr.puts "After culling for uniqueness, #{vector_compatible_cuts.length} vector-compatible sequences remain"
  # Remove ones to specifically avoid because they are expensive or unavailable or whatever
  vector_compatible_cuts = vector_compatible_cuts.reject{|c| options[:enzymes_to_avoid].include?(c.enzyme)}
  $stderr.puts "After culling those enzymes that are to be globally avoided, #{vector_compatible_cuts.length} vector-compatible sequences remain"

  # ============================================================================
  # Attempt to choose a suitable enzyme
  # first, has the enzyme to be used been specified by the user already?
  enzyme = nil
  if options[:enzyme]
    hits = vector_compatible_cuts.select do |v|
      v.enzyme == options[:enzyme]
    end
    if hits.length > 1
      raise Exception, "what the hell? This shouldn't happen. Too many enzymes found matching"
    elsif hits.length == 0
      $stderr.puts "You specified to use enzyme `#{options[:enzyme]}', but that enzyme wasn't found to be vector-compatible. Choose another?"
      exit
    else
    # all is well.
      enzyme = hits[0]
      $stderr.puts "Using enzyme #{enzyme} as requested"
    end
  else
  # if there's one in the freezer, the choice is clear even though we haven't specified things manually
    if freezer_compatible_cuts.length > 0
      # choose the one with the least length for easier cloning
      enzyme = freezer_compatible_cuts.max{|a,b| a.stop <=> b.stop}
      $stderr.puts "Found a compatible enzyme in the freezer: #{enzyme.to_s} (All freezer-compatible choices: #{freezer_compatible_cuts.collect{|e| e.to_s}.join(', ')} but I've chosen the one with the shortest amount of sequence)"
    else
    # if there are no freezer-compatible cuts, then need to ask the user for which enzyme we are going to use, and that means re-running the whole script
      $stderr.puts "No compatible enzymes found, so I need your help. Tell me the name of the enzyme that you want me to use (it must be in the vector compatible list), using the -e option to this script."
      $stderr.puts "Here's a list:"
      vector_compatible_cuts.each do |v|
        $stderr.puts v
      end
    end
  end
  Process.exit if enzyme.nil? #don't really understand why this is required, but
  # otherwise errors happen if no compatible enzymes are found
  $stderr.puts "Using enzyme #{enzyme.enzyme}"

  # ============================================================================
  # Attempt to design primers

  # If the last 3 bases translate to a stop codon, then remove them from the
  # sequence
  if Bio::Sequence::NA.new(upstream_sequence[upstream_sequence.length-3..upstream_sequence.length-1]).translate == '*'
    $stderr.puts "Found a stop codon in the last position of the sequence. Removing this."
  upstream_sequence = upstream_sequence[0..upstream_sequence.length-4]
  end

  $stderr.puts "Attemping to find primers in the following sequence, with a cut site at #{enzyme.start} with #{enzyme.enzyme}"
  #$stderr.puts upstream_sequence

  # Use primer 3 to find a primer under various constraints. Start with the
  # default GC clamp, and then keep going down to options[:real_gc_minimum]
  current_gc_clamp = options[:primer_gc_clamp]
  found_primer = false
  while !found_primer and current_gc_clamp >= options[:real_gc_clamp_minimum]
    Tempfile.open('primer3input') do |tempfile|
      min_length = upstream_sequence.length-enzyme.start+options[:five_prime_homologous_recombination_minimum]
      max_length = upstream_sequence.length-enzyme.start+options[:five_prime_homologous_recombination_minimum]+options[:primer_search_area]
      input_hash = {
        'SEQUENCE_TEMPLATE' => upstream_sequence,
        # needs to be a certain size - needs to be a certain amount 5' of the
        # enzyme cut site for smooth homologous recombination
        'PRIMER_PRODUCT_SIZE_RANGE' => "#{min_length}-#{max_length}",
        'PRIMER_GC_CLAMP' => current_gc_clamp,
        'PRIMER_PAIR_MAX_DIFF_TM'=>1,
        # 'PRIMER_TM_SANTALUCIA' => 1, #as recommended by primer3 - seems not to
        # work for some reason.
        'PRIMER_SALT_CORRECTIONS' => 1, #as recommended by primer3
        'SEQUENCE_FORCE_RIGHT_START' => upstream_sequence.length-1, #right primer
        # must cover the last base, which is the base before the stop codon
        'PRIMER_EXPLAIN_FLAG'=>1, #be verbose
        'PRIMER_MIN_SIZE'=>15,
        'PRIMER_MAX_SIZE'=>35,
        'PRIMER_MIN_TM'=>50.0,
        'PRIMER_MAX_TM'=>75.0,
        #'PRIMER_LEFT_MAX_POLY_X' => 3, #left primer shouldn't have runs in it
        # #doesn't seem to accept this, for some reason
        #'PRIMER_LEFT_MAX_SELF_ANY' => 3, #doesn't seem to accept this, for some
        # reason
        'PRIMER_PAIR_MAX_COMPL_ANY' => 3, #primer dimers
        # 'PRIMER_THERMODYNAMIC_ALIGNMENT' => 1, # commented out because it didn't seem to give good primers
        # 'PRIMER_THERMODYNAMIC_PARAMETERS_PATH' => '/home/ben/bioinfo/primer3-2.2.3/src/primer3_config/',
      }
      record = BoulderIO::Record.new(input_hash)
      tempfile.print record.to_s
      tempfile.close

      puts "Querying primer3 with the following input: "
      print record.to_s

      #run primer3
      Tempfile.open('primer3out') do |primer3out|
        puts `primer3_core -strict #{tempfile.path} >#{primer3out.path}`
        #puts `cat #{primer3out.path}`

        # Was there a primer found?
        # Read the output file, which is in boulder I/O format
        result = Primer3Result.create_from_primer3_output_filename(primer3out.path)

        if result.yeh?
        	puts "Using enzyme #{enzyme.enzyme}, which cuts at #{enzyme.start}"
          (0..4).each do |primer_number|
            puts "Found primer #{primer_number}: left: #{result["PRIMER_LEFT_#{primer_number}_SEQUENCE"]}"
            puts "Found primer #{primer_number}: right: #{result["PRIMER_RIGHT_#{primer_number}_SEQUENCE"]}"
            puts "Found primer #{primer_number}: product size: #{result["PRIMER_PAIR_#{primer_number}_PRODUCT_SIZE"]}"
          end
        found_primer = true #success!
        else
          puts "No suitable primers. Some possible reasons:"
          puts "PRIMER_ERROR: #{result['PRIMER_ERROR']}"
          puts "PRIMER_PAIR_EXPLAIN: #{result['PRIMER_PAIR_EXPLAIN']}"
          puts "PRIMER_RIGHT_EXPLAIN: #{result['PRIMER_RIGHT_EXPLAIN']}"
          puts "PRIMER_LEFT_EXPLAIN: #{result['PRIMER_LEFT_EXPLAIN']}"
        current_gc_clamp -= 1 #relax the GC clamp constraint
        end
      end
    end
  end
end
