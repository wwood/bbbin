#!/usr/bin/env ruby

# Takes in a reference fasta file and a file of SNPs from dnadiff output e.g. out.snps

# Outputs information about the snps

require 'bio'
require 'csv'
require 'optparse'
require 'bio-samtools'
require 'bio-logger'
$:.unshift File.join([ENV['HOME'], %w(git bioruby-pileup_iterator lib)].flatten)
require 'bio-pileup_iterator'
require 'pp'

class SNP
  attr_reader :ref, :aligned
  attr_accessor :position, :ref_name, :contig_name, :contig_position, :opposing_direction
  
  def ref=(ref)
    if ref == '.'
      @ref = nil
    else
      @ref = ref
    end
  end
  
  def aligned=(aligned)
    if aligned == '.'
      @aligned = nil
    else
      @aligned = aligned
    end
  end
  
  def type
    return 'insertion' if @ref.nil? and !@aligned.nil?
    return 'deletion' if deletion?
    return 'substitution' if !@ref.nil? and !@aligned.nil?
    p inspect
    raise
  end
  
  def deletion?
    !@ref.nil? and @aligned.nil?
  end
  
  def homopolymer_length(seq)
    count = 1 #count itself
    
    # How many out in front?
    i = 1
    next_base = seq[@position+i]
    while next_base == ref
      count += 1
      i += 1
      next_base = seq[@position+i]
    end

    # How many out back?
    i = -1
    next_base = seq[@position+i]
    while next_base == ref
      count += 1
      i -= 1
      next_base = seq[@position+i]
    end
    
    return count
  end
end




class Bio::DB::Sam
      #calls the mpileup function, opts is a hash of options identical to the command line options for mpileup.
      #is an iterator that yields a Pileup object for each postion
      #the command line options that generate/affect BCF/VCF are ignored ie (g,u,e,h,I,L,o,p)
      #call the option as a symbol of the flag, eg -r for region is called :r => "some SAM compatible region"
      #eg bam.mpileup(:r => "chr1:1000-2000", :q => 50) gets the bases with quality > 50 on chr1 between 1000-5000
      def mpileup( opts )

              raise SAMException.new(), "No BAMFile provided" unless @sam and @binary
              raise SAMException.new(), "No FastA provided" unless @fasta_path
              #long option form to short samtools form..
              long_opts = {
              :region => :r,
              :illumina_quals => :six,
              :count_anomalous => :A,
              :no_baq => :B,
              :adjust_mapq => :C,
              :max_per_bam_depth => :d,
              :extended_baq => :E,
              :exclude_reads_file => :G,
              :list_of_positions => :l,
              :mapping_quality_cap => :M,
              :ignore_rg => :R,
              :min_mapping_quality => :q,
              :min_base_quality => :Q
              }
              ##convert any long_opts to short opts
              temp_opts = opts.dup
              opts.each_pair do |k,v|
                if long_opts[k]
                  temp_opts[long_opts[k]] = v
                  temp_opts.delete(k)
                end
              end
              opts = temp_opts
              ##remove any calls to -g or -u for mpileup, bcf output is not yet supported
              ##and also associated output options
              [:g, :u, :e, :h, :I, :L, :o, :p].each {|x| opts.delete(x) }

              sam_opts = []
              #strptrs << FFI::MemoryPointer.from_string("mpileup")
              opts.each do |k,v|
                next unless opts[k] ##dont bother unless the values provided are true..
                k = '6' if k == :six
                k = '-' + k.to_s
                sam_opts << k #strptrs << FFI::MemoryPointer.from_string(k)
                sam_opts << v.to_s unless ["-R", "-B", "-E", "-6", "-A"].include?(k) #these are just flags so don't pass a value... strptrs << FFI::MemoryPointer.from_string(v.to_s)
              end
              sam_opts = sam_opts + ['-f', @fasta_path, @sam]
              sam_command = "samtools mpileup #{sam_opts.join(' ')} 2>/dev/null"

              sam_pipe = IO.popen(sam_command)
              lines = sam_pipe.readlines
              lines.each do |line|
                yield Bio::DB::Pileup.new(line)
              end
              sam_pipe.close
              return lines
              #strptrs << FFI::MemoryPointer.from_string('-f')
              #strptrs << FFI::MemoryPointer.from_string(@fasta_path)
              #strptrs << FFI::MemoryPointer.from_string(@sam)
              #strptrs << nil

              # Now load all the pointers into a native memory block
              #argv = FFI::MemoryPointer.new(:pointer, strptrs.length)
              #strptrs.each_with_index do |p, i|
              # argv[i].put_pointer(0, p)
              #end

              #old_stdout = STDOUT.clone
              #read_pipe, write_pipe = IO.pipe()
              #STDOUT.reopen(write_pipe)
                #int bam_mpileup(int argc, char *argv[])
               # Bio::DB::SAM::Tools.bam_mpileup(strptrs.length - 1,argv)
                #if fork
                # write_pipe.close
                # STDOUT.reopen(old_stdout) #beware .. stdout from other processes eg tests calling this method can get mixed in...
                # begin
                # while line = read_pipe.readline
                # yield Pileup.new(line)
                # end
                # rescue EOFError
                # read_pipe.close
                # Process.wait
                # end
                #end
            end
end






SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')
# Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
Bio::Log::CLI.logger('stderr'); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)


# Parse options
OVERALL = 'overall'
DELETION = 'deletion'
DELETION_STRANDEDNESS = 'deletion_strandedness'

options = {
  :operation => OVERALL,
}
o = OptionParser.new do |opts|
  opts.banner = [
    'Usage: classify_iontorrent_snps.rb <reference fasta file> <dnadiff out.snps>',
    "\tTakes in a reference fasta file and a file of SNPs from dnadiff output e.g. out.snps\n"+
    "\tOutputs information about the snps\n"
  ]
  opts.on("-o", "--overall", "Overall classification of SNPs IN/DEL/SUB breakdown") do |f|
    options[:operation] = OVERALL
  end
  opts.on('-d', "--deletion", "Investigate the deletions specifically") do |e|
    options[:operation] = DELETION
  end
  opts.on("--deletion-strandedness", "Investigate the strand bias at deletion positions specifically. Requires --pileup") do |e|
    options[:operation] = DELETION_STRANDEDNESS
  end
  # opts.on("--ace ACE_FILE", "ACE file representing the assembly") do |e|
    # options[:ace_file] = e
  # end
  opts.on("--bam BAMFILE", "(indexed) BAM file representing the assembly") do |e|
    options[:bam_file] = e
  end
end
o.parse!
if ARGV.length != 2
  $stderr.puts o.banner
  exit
end



# Read in the SNPs
log.info "Reading SNPs"
snps = []
CSV.foreach(ARGV[1],:col_sep => "\t") do |row|
  snp = SNP.new
  
  # [P1]  [SUB] [SUB] [P2]  [BUFF]  [DIST]  [LEN R] [LEN Q] [FRM] [TAGS]
  snp.ref = row[1]
  snp.aligned = row[2]
  snp.position = row[0].to_i
  snp.contig_position = row[3].to_i
  if row[9] == '-1'
    snp.opposing_direction = true
  else
    snp.opposing_direction = false
  end
  snp.ref_name = row[10]
  snp.contig_name = row[11]
  snps.push snp
end

types = {}
snps.each do |snp|
  types[snp.type] ||= 0
  types[snp.type] += 1
end

log.info "Read in #{snps.length} SNPs"


# Classify each of SNPs
if options[:operation] == OVERALL
  types.each do |type, count|
    puts [type,count].join("\t")
  end



# Print out more specifics on the deletions
elsif options[:operation] == DELETION
  # Read in the fasta sequences
  seqs = {}
  Bio::FlatFile.foreach(ARGV[0]){|s| seqs[s.definition.split(' ')[0]] = s}
  
  #Print out the surrounding sequences around the deletions
  snps.each do |snp|
    next unless snp.deletion?
    
    ref = seqs[snp.ref_name]
    surround4 = ref.seq[snp.position-2..snp.position+1]
    homopolymer_length = snp.homopolymer_length(ref.seq) 
    
    puts [
      snp.ref_name,
      snp.position,
      snp.ref,
      surround4,
      homopolymer_length,
    ].join("\t")
  end
  
  
elsif options[:operation] = DELETION_STRANDEDNESS
  # Read in the fasta sequences
  seqs = {}
  Bio::FlatFile.foreach(ARGV[0]){|s| seqs[s.definition.split(' ')[0]] = s}
  
  gag_type_contexts = %w(GAAG CTTC  AGGC GCCT  GCCG CGGC  GCCA TGGC)
  
  # headers
  puts %w(
  ref_name
  contig_name
  ref_position
  contig_position
  ref_base
  context
  forwards_insertion
  forwards_no_insertion
  backwards_insertion
  backwards_no_insertion
  opposing_direction
  ).join("\t")
  
  #Print out the surrounding sequences around the deletions
  snps.each do |snp|
    next unless snp.deletion?
    
    ref = seqs[snp.ref_name]
    surround4 = ref.seq[snp.position-2..snp.position+1]
    
    # Disregard non-gag types
    next unless gag_type_contexts.include?(surround4)
    
    # How many reads are 1) in each direction?
    # get pileup for this range using samtools
    bam = Bio::DB::Sam.new(:bam => options[:bam_file], :fasta => ARGV[0])
    region = nil
    if snp.opposing_direction
      region = "'#{snp.contig_name}:#{snp.contig_position}-#{snp.contig_position+1}'"
    else
      region = "'#{snp.contig_name}:#{snp.contig_position+1}-#{snp.contig_position+2}'"
    end
    
    pileups = bam.mpileup(:r => region){}
    iterator = Bio::DB::PileupIterator.new(StringIO.new(pileups[0]))
    forwards_insertion = 0
    forwards_no_insertion = 0
    backwards_insertion = 0
    backwards_no_insertion = 0
    
    iterator.each do |pileup|
      pileup.reads.each do |read|
        has_insertion = read.insertions.length > 0
        
        if read.direction == '+' and has_insertion
          forwards_insertion += 1
        elsif read.direction == '+' and !has_insertion
          forwards_no_insertion += 1
        elsif read.direction == '-' and has_insertion
          backwards_insertion += 1
        elsif read.direction == '-' and !has_insertion
          backwards_no_insertion += 1
        else
          log.warn "Unexpected orientation found in pileup #{region} '#{read.direction}', ignoring"
        end
      end
      break #only want the first position
    end
    
    
    puts [
      snp.ref_name,
      snp.contig_name,
      snp.position,
      snp.contig_position,
      snp.ref,
      surround4,
      forwards_insertion,
      forwards_no_insertion,
      backwards_insertion,
      backwards_no_insertion,
      snp.opposing_direction,
    ].join("\t")
  end
end






