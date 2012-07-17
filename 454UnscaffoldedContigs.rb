#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'csv'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

module Bio
  class Newbler
    class Contig
      attr_accessor :length, :coverage
    end
    
    
    # Position Consensus Quality_Score Depth Signal  StdDeviation
    # >contig00001  1
    # 1 C 64  6 2.04  0.60
    # 2 C 20  6 2.04  0.60
    # 3 G 64  6 1.96  0.55
    # 4 G 20  6 1.96  0.55
    class AlignmentInfoFile
      # Iterate through each row, returning contig_name, columns (string then Array of strings)
      def self.foreach(filename)
        current_contig = nil
        
        CSV.foreach(filename, :header => true, :col_sep => "\t") do |row|
          if matches = row[0].match(/^>(contig\d{5})$/)
            current_contig = matches[1]
          else
            yield current_contig, row
          end
        end
      end
      
      # Return a hash of contig name to Contig objects with filled in coverage and length details
      def self.hash(filename)
        hash = {}
        add_contig = lambda do |contig_name, total_count, length|
          contig = Contig.new
          contig.length = length
          contig.coverage = total_count.to_f/length
          hash[contig_name] = contig
        end
        
        foreach(filename) do |contig, cols|
          last_contig = contig if last_contig.nil?
          
          if contig==last_contig 
            count += cols[3].to_i
            total_length += 1
          else
            add_contig.call(last_contig, count, total_length)
          end
        end
        add_contig.call(last_contig, count, total_length)
        
        return hash
      end
    end
  end
end

if __FILE__ == $0 #needs to be removed if this script is distributed as part of a rubygem
  # Parse command line options into the options hash
  options = {
    :logger => 'stderr',
  }
  o = OptionParser.new do |opts|
    opts.banner = "
      Usage: #{SCRIPT_NAME}
      
      Returns a list of contigs with some minimum coverage, and that are not scaffolded at all by newbler.\n\n"
    # Example option
    opts.on("-s", "--scaffolds-txt 454SCAFFOLDS.TXT", "From newbler [required]") do |f|
      options[:scaffolds_file] = f
    end
    opts.on("-m", "--min-coverage MIN_COVERAGE", "From newbler [required]") do |f|
      options[:min_coverage] = f.to_f
      raise if options[:min_coverage] < 0 
    end
    opts.on("-c", "--contig-graph 454CONTIGGRAPH.TXT", "From newbler [required unless -a/--alignment-info is defined]") do |f|
      options[:contigs_graph_file] = f
    end
    opts.on("-a", "--alignment-info 454ALIGNMENT_INFO.TXT", "From newbler [required unless -c/--contig-graph is defined]") do |f|
      options[:alignment_info_file] = f
    end
    
    # logger options
    opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") do |q|
      Bio::Log::CLI.trace('error')
    end
    opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") do | name |
      options[:logger] = name
    end
    opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG") do | s |
      Bio::Log::CLI.trace(s)
    end
  end
  o.parse!
  if ARGV.length != 0
    $stderr.puts o
    exit 1
  end
  # Setup logging
  Bio::Log::CLI.logger(options[:logger]) #bio-logger defaults to STDERR not STDOUT, I disagree
  log = Bio::Log::LoggerPlus.new(LOG_NAME)
  Bio::Log::CLI.configure(LOG_NAME)
  
  
  # Read in the scaffolds file, these are the ones that should not be printed out
  scaffolded_contigs = {}
  CSV.foreach(options[:scaffolds_file], :col_sep => "\t") do |row|
    scaffolded_contigs[row[5]] = true if row[4]=='W'
  end
  log.info "Found #{scaffolded_contigs.length} contigs in scaffolds"
  
  low_coverage_count = 0
  ignored_contigs_count = 0
  printed_contigs_count = 0
  
  process_contig = lambda do |contig_name, coverage|
    if coverage < options[:min_coverage]
      low_coverage_count += 1
    elsif scaffolded_contigs[contig_name]
      ignored_contigs_count += 1
    else
      printed_contigs_count += 1
      puts contig_name
    end
  end
  
  # Either use the contigs graph file or the alignment info file
  if options[:contigs_graph_file]
    CSV.foreach(options[:contigs_graph_file], :col_sep => "\t") do |row|
      break if row[0] == 'C'
      raise unless row[0].to_i > 0 #indicates unexpected file structure
      
      coverage = row[3].to_f
      contig_name = row[1]
      process_contig.call(contig_name, coverage)
    end
    
  else
    last_contig = nil
    count = 0
    total_length = 0
    hash = Bio::Newbler::AlignmentInfoFile.hash(options[:alignment_info_file])
    hash.each do |name, contig|
      process_contig.call(name, contig.coverage)
    end
  end
  
  log.info "Ignored due to low coverage: #{low_coverage_count}"
  log.info "Ignored due to being in a scaffold: #{ignored_contigs_count}"
  log.info "Printed contigs: #{printed_contigs_count}"

end















