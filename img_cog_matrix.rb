#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'csv'
require 'progressbar'

$:.push File.join(ENV['HOME'],'git','bioruby-taxonomy_definition_files','lib')
require 'bio-taxonomy_definition_files' #has IMG taxonomy parser file

if __FILE__ == $0 #needs to be removed if this script is distributed as part of a rubygem
  SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')
  
  # Parse command line options into the options hash
  options = {
    :logger => 'stderr',
    :img_base_directory => '/srv/whitlam/bio/db/img/4.0/genomes/all',
    :img_metadata_file => '/srv/whitlam/bio/db/img/4.0/metadata/img_metadata_4_0_FIXED.csv',
  }
  o = OptionParser.new do |opts|
    opts.banner = "
      Usage: #{SCRIPT_NAME}
      
      Output a matrix of IMG COGs for each species, as well as their temperature range, and domain\n\n"
      

    # logger options
    opts.separator "\n\tVerbosity:\n\n"
    opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {Bio::Log::CLI.trace('error')}
    opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
    opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| Bio::Log::CLI.trace(s)}
  end; o.parse!
  if ARGV.length != 0
    $stderr.puts o
    exit 1
  end
  # Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
  Bio::Log::CLI.logger(options[:logger]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)
  
  
  # Read in the taxonomy file
  taxonomies = Bio::IMG::TaxonomyDefinitionFile.read(options[:img_metadata_file])
  raise unless taxonomies.length > 1000
  log.info "Read in #{taxonomies.length} different taxonomy entries"
  
  taxon_cogs = {}
  progress = ProgressBar.new('cog_caching',taxonomies.length)
  
  taxonomies.each do |taxonomy|
    cog_file = File.join options[:img_base_directory], taxonomy.taxon_id.to_s, "#{taxonomy.taxon_id}.cog.tab.txt"
    if File.exist?(cog_file)
      cogs = {}
      CSV.foreach(cog_file, :col_sep => "\t", :headers => true) do |row|
        cog = row[9]
        cogs[cog] ||= 0
        cogs[cog] += 1
      end
      
      raise "Duplicate taxon id" if taxon_cogs[taxonomy]
      taxon_cogs[taxonomy] = cogs
    else
      log.warn "Didn't find any COG file for #{taxonomy.taxon_id} #{taxonomy.genus_species}"
    end
    progress.inc
  end
  progress.finish
  
  # Print results. Get a list of all COGs
  all_cogs = taxon_cogs.collect{|taxonomy, cog_counts| cog_counts.keys}.flatten.uniq.sort
  puts [
         'taxon_id',
         'genus',
         'species',
         'temperature_range',
         all_cogs].flatten.join "\t"
  taxon_cogs.each do |taxonomy, cog_counts|
    print [
      taxonomy.taxon_id,
      taxonomy.genus,
      taxonomy.species,
      taxonomy.attributes['Temperature Range'],
    ].join "\t"
    
    all_cogs.each do |cog_id|
      print "\t"
      if cog_counts[cog_id].nil?
        print 0
      else
        print cog_counts[cog_id]
      end
    end
    puts
  end
end #end if running as a script
