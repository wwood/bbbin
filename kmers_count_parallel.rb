#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'tmpdir'
require 'tempfile'
require 'bio-commandeer'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',

  :chunk_size => 5000000,
  :num_threads => 44,
  :threshold => 3,
  :kmer => 51,
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <arguments>

    Description of what this program does...\n\n"

  opts.on("-i", "--input FILE",String,"Input fastq.gz file [required]") do |arg|
    options[:input_file] = arg
  end
  opts.on("-o", "--output FILE",String,"Output space separated kmer, fastq.gz file [required]") do |arg|
    options[:output_file] = File.absolute_path arg
  end

  opts.separator "\nOther optional parameters:\n\n"
  opts.on("-k", "--kmer INTEGER",Integer,"kmer size to count [default: #{options[:kmer] }") do |arg|
    options[:kmer] = arg
  end
  opts.on("--threads INTEGER",Integer,"Number of threads. Note the number of CPUs used will vary because of piping. [default: #{options[:num_threads] }") do |arg|
    options[:num_threads] = arg
  end
  opts.on("--chunk-size INTEGER",Integer,"Number of lines in fastq file that each chunk is. Must be divisible by 4. [default: #{options[:chunk_size] }") do |arg|
    options[:chunk_size] = arg
  end
  opts.on("--threshold INTEGER",Integer,"Ignore kmers with less or equal abundance count than this, so the output file is reduced in size [default: #{options[:threshold] }") do |arg|
    options[:threshold] = arg
  end

  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end; o.parse!
if ARGV.length != 0 or options[:input_file].nil? or options[:output_file].nil?
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME); log.outputters[0].formatter = Log4r::PatternFormatter.new(:pattern => "%5l %c %d: %m", :date_pattern => '%d/%m %T')

input_path_absolute = File.absolute_path options[:input_file]

Dir.mktmpdir do |tmpdir|
  #tmpdir = '/tmp/test'
  Dir.chdir(tmpdir) do
    # Copy reads file to temporary directory
    log.info "Copying raw reads file `#{input_path_absolute}' into temporary directory.."
    local_input_name = 'input.fq.gz'
    FileUtils.cp input_path_absolute, local_input_name
    log.info "Finished copying reads file."

    #log.info "Counting the number of lines in the file..."
    #num_input_lines = Bio::Commandeer.run "pigz -cd '#{local_input_name }' |wc -l", :log => log
    #num_input_lines = num_input_lines.strip.to_i
    #log.info "Found #{num_input_lines } lines in the input file"

    # Split it into chunks
    Dir.mkdir 'chunks'
    #chunk_size = options[:chunk_size]
    #num_chunks = num_input_lines / chunk_size +1
    log.info "Splitting the reads file up into chunks.."
    cmd = "pigz -cd #{local_input_name} |split -l #{options[:chunk_size] } --filter 'pigz > $FILE.gz' --additional-suffix .fq - chunks/chunk_"
    Bio::Commandeer.run cmd, :log => log
    log.info "Finished chunking"

    # each chunk in parallel: decompress, kmers_count, sort, recompress
    chunk_files = Dir.glob('chunks/*')
    sort_buffer_size = 100 / options[:num_threads]
    sort_buffer_size = 1 if sort_buffer_size == 0
    kmer_count_cmds = chunk_files.collect do |chunkf|
      "<(kmers_count -r -q -k #{options[:kmer] } <(pigz -cd #{chunkf} ) |sort -k 1b,1 --parallel=1 --buffer-size='#{sort_buffer_size}%')"
    end

    # join them all up together, remove low coverage kmers, output the file compressed into the output file
    join_args = 'join -a1 -a2 -e0'
    log.info "kmer_counting, joining, awk'ing and recompressing #{kmer_count_cmds.length} chunks.."
    log.debug "chunk_files: #{kmer_count_cmds.inspect}"
    #     join -a1 -a2 -e 0 -o0,1.2,2.2 <(pigz -cd chunk_counts/20110816_S1D.1.chunk1.fq.gz.counts.gz) <(pigz -cd chunk_counts/20110816_S1D.1.chunk5000001.fq.gz.counts.gz) |\
    # join -a1 -a2 -e 0 -o0,1.`seq 2 3 |tr '\n' : |sed 's/:/,1./g' |sed 's/...$//'`,2.2 - <(pigz -cd chunk_counts/20110816_S1D.1.chunk10000001.fq.gz.counts.gz) |\
    # join -a1 -a2 -e 0 -o0,1.`seq 2 4 |tr '\n' : |sed 's/:/,1./g' |sed 's/...$//'`,2.2 - <(pigz -cd chunk_counts/20110816_S1D.1.chunk15000001.fq.gz.counts.gz) |\
    cmd = "#{join_args} -o0,1.2,2.2 #{kmer_count_cmds[0] } #{kmer_count_cmds[1] }"
    # add all the join lines
    (2...kmer_count_cmds.length).each do |i|
      outputs = (2..(i+1)).collect{|j| "1.#{j}"}.join(',')
      cmd += " | #{join_args} -o0,#{outputs},2.2 - #{kmer_count_cmds[i] }"
    end
    # add all the stuff afterwards
    # |awk "(\$`seq 2 20 |tr '\n' ':' |sed 's/.$//' |sed 's/:/+\$/g'`>3){print \$1,\$`seq 2 20 |tr '\n' ':' |sed 's/.$//' |sed 's/:/+\$/g'`}" |pigz >20110816_S1D.1.fq.kmer51.gz
    sum = (2..(kmer_count_cmds.length+1)).collect{|j| "\\\$#{j}"}.join('+')
    cmd += " | parallel --gnu --pipe --keep-order --block 100M awk \"'(#{sum}>#{options[:threshold] }){print \\$1,#{sum}}'\" |pigz >#{options[:output_file] }"
    Tempfile.open('bash_script') do |tf|
    #File.open('bash_script','w') do |tf|
      tf.puts cmd
      tf.close
      log.info "Running script #{tf.path} with contents #{File.open(tf).read.inspect}"

      Bio::Commandeer.run "bash #{tf.path}", :log => log
    end
    log.info "Wrote output to #{options[:output_file] }"
  end
end


