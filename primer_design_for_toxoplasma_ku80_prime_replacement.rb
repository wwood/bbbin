#!/usr/bin/env ruby

# A script to make the design primers for T. gondii ku80 3' replacement. Currently
# this script only works by taking a sequence, and it spits out a list of
# possible enzyme digestion sites. This can be helpful in itself. It then
# tries to find primer sites, but I wouldn't trust those results just yet, and
# it isn't likely to find any suitable primers anyway.
#
# The list of enzymes that are compatible can be changed by modifying this script,
# or using the -e option. By default, it is overly ben-specific. 

require 'rubygems'
require 'bio'
require 'optparse'
require 'eupathdb_gene_information_table' #from the reubypathdb rubygem

require 'tempfile'
require '/home/ben/bin/primer_design_extras' # a loose collection of helpful

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
    :print_vector_compatible_sites => false,
    :print_freezer_compatible_sites => false,
#    :recombination_sequence_buffer => 350, #how far way from priming sites can
    # the restriction site be?
    :three_prime_homologous_recombination_minimum => 800,
    :five_prime_homologous_recombination_minimum => 350,
    :enzyme => nil,
    :primer_temperature_minimum => 53,
    :primer_temperature_optimal => 60,
    :primer_temperature_maximum => 75,
    :primer_gc_clamp => 2,
    :enzymes_on_order => [], #extra enzymes that will be in the freezer soon
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
  # Find the unique restriction sites that are in that 1kb, and a far enough
  # distance away from the 3' end
  Tempfile.open('first_fasta') do |first_fasta_file|
    first_fasta_file.puts '>seq'
    first_fasta_file.puts upstream_sequence
    first_fasta_file.close
    Tempfile.open('restrict') do |restrict_output_file|
    # find the unique restriction sites
      cmd = "restrict -sequence #{first_fasta_file.path} -sitelen 4 -single -enzymes all -outfile #{restrict_output_file.path}"
      puts `#{cmd}`

      # extract the list of pure sites
      Tempfile.open('awk_output') do |awk_out|
        cmd = "awk 'NF==9{print $0}' #{restrict_output_file.path} >#{awk_out.path}"
        `#{cmd}`
        print 'Found this many unique sites in the 1500bp: '
        puts `wc -l #{awk_out.path} |awk '{print $1}'`

        # merge the list of enzymes with the ones that are known to not be
        # a) commercially available
        # b) don't also cut the LIC Ku80 vector
        Tempfile.open('enzyme_match1') do |enzyme_match1|
          `enzyme_matcher.rb -f '#{options[:possible_enzymes_to_cut_with_filename]}' #{awk_out.path} >#{enzyme_match1.path}`
          print 'found this many vector-compatible: '
          puts `wc -l #{enzyme_match1.path} |awk '{print $1}'`
          puts `cat #{enzyme_match1.path}` if options[:print_vector_compatible_sites]

          # parse restriction sites
          vector_compatible_cuts = enzyme_cuts_array_from_filename.call(enzyme_match1.path)

          # Do we have these enzymes on the lab's list of enzymes in the freezer?
          Tempfile.open('enzyme_match2') do |enzyme_match2|
            args = ''
            unless options[:enzymes_on_order].nil?
              args += "-e #{options[:enzymes_on_order].join(',')}"
            end
            `enzyme_matcher.rb -f '#{options[:enzymes_in_freezer_filename]}' #{args} #{enzyme_match1.path} >#{enzyme_match2.path}`
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

  # Are there any freezer-enzymes that are far enough away from the 3' end?
  freezer_compatible_cuts = freezer_compatible_cuts.select do |f|
    f.stop < upstream_sequence.length-options[:three_prime_homologous_recombination_minimum]
  end
  enzyme = nil

  # first, has the enzyme to be used been specified by the user already?
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

  $stderr.puts "Attemping to find primers in the following sequence:"
  $stderr.puts upstream_sequence

  # First step, find a right primer that is just less than 72 degrees
  o = OligoDesigner.new
  revcomp = Bio::Sequence::NA.new(upstream_sequence).reverse_complement.to_s

  # Get a list of oligos that are possible given the constraints
  possible_right_primers = o.order(revcomp,
  options[:primer_temperature_minimum],
  options[:primer_temperature_optimal],
  options[:primer_temperature_maximum],
  options[:primer_gc_clamp]
  )
  $stderr.puts "Found #{possible_right_primers.length} right primers within constraints. Now attempting to design primers with these.."

  first = true

  possible_right_primers.each do |right_primer_sequence|

    $stderr.puts "Trying right primer #{right_primer_sequence}, melting temp #{o.melting_temperature right_primer_sequence}."
    Tempfile.open('primer3input') do |tempfile|
      input_hash = {
        'PRIMER_RIGHT_INPUT' => right_primer_sequence,
        'SEQUENCE' => upstream_sequence,
        # needs to be a certain size - needs to be a certain amount 5' of the enzyme cut site for smooth homologous recombination
        'PRIMER_PRODUCT_SIZE_RANGE' => "#{upstream_sequence.length-enzyme.start+options[:five_prime_homologous_recombination_minimum]}-#{upstream_sequence.length}",
        'PRIMER_GC_CLAMP' => options[:primer_gc_clamp],
      }
      record = BoulderIO::Record.new(input_hash)
      tempfile.print record.to_s
      tempfile.close

      if first
        $stderr.puts "Querying primer3 with the following input: "
      $stderr.puts record.to_s
      first = false
      end

      #run primer3
      Tempfile.open('primer3out') do |primer3out|
        puts `primer3_core -strict_tags <#{tempfile.path} >#{primer3out.path}`

        # Was there a primer found?
        # Read the output file, which is in boulder I/O format
        result = Primer3Result.create_from_primer3_output_filename(primer3out.path)

        if result.yeh?
          $stderr.puts "Found primer: #{result['PRIMER_LEFT_SEQUENCE']}"
        # print out the other primers if more than one was found.
        # TODO
        else
          $stderr.puts "No suitable primers"
        end
      end
    end
  # right primer: using the reverse complement of the first however many bases
  # left primer: try to design this
  # ensure that the length of the product is at least the distance from the  end of the sequence to the restriction site + 350bp

  # If primers found, output, otherwise
  # a) try another restriction enzyme
  # b) try with more than 1kb being extracted
  end

end