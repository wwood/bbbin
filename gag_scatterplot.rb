#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'forkmanager'
require 'bio'
require 'bio-samtools'
require 'progressbar'

$:.unshift File.join([ENV['HOME'], %w(git bioruby-pileup_iterator lib)].flatten)
require 'bio-pileup_iterator'


module Bio::Sequence::Common
  def window_search(window_size, step_size = 1)
    last_step = 0
    0.step(self.length - window_size, step_size) do |i| 
      yield self[i, window_size], i+1                        
      last_step = i
    end                          
    return self[last_step + window_size .. -1] 
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


if __FILE__ == $0 #needs to be removed if this script is distributed as part of a rubygem
  SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')
  
  # Parse command line options into the options hash
  options = {
    :logger => 'stderr',
    :contexts_to_interrogate => %w(GAAG CTTC  AGGC GCCT  GCCG CGGC  GCCA TGGC),
    # :num_threads => 10,
  }
  o = OptionParser.new do |opts|
    opts.banner = "
      Usage: #{SCRIPT_NAME} <arguments>
      
      Take a genome, an assembly, and a list of 4mers. Output for each position that matches the 4mer (as a gag error) and is mapped to part of the assembly, the number of reads at that part that have deletions, and those that don\'t\n\n"
      
    opts.on("--bam BAMFILE_NAME", "samfile of the reads mapped to the reference [required]") do |arg|
      options[:bam_file] = arg
    end
    opts.on("-r", "--reference FASTA_FILE", "Fasta file of the reference assembly [required]") do |arg|
      options[:reference_fasta_file] = arg
    end
    # opts.on("-t", "--threads NUM_THREADS", "Threads to work with [default: #{options[:num_threads]}") do |arg|
      # options[:num_threads] = arg.to_i
    # end
    
    # logger options
    opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {Bio::Log::CLI.trace('error')}
    opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
    opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| Bio::Log::CLI.trace(s)}
  end
  o.parse!
  if ARGV.length != 0
    $stderr.puts o
    exit 1
  end
  # Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
  Bio::Log::CLI.logger(options[:logger]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)
  
   
  sam = Bio::DB::Sam.new(
    :fasta => options[:reference_fasta_file],
    :bam => options[:bam_file]
  )
  
  progress = nil
  # data structure retrieval and handling
  # pm.run_on_finish { # called BEFORE the first call to start()
      # |pid,exit_code,ident,exit_signal,core_dump,data_structure|
     # # retrieve data structure from child
     # #puts data_structure.inspect
     # puts data_structure unless data_structure=={}
     # progress.inc
  # }
  # Print headers
  puts %w(
    reference
    position
            seq
            coverage
            num_reads_with_deletion
            num_reads_without_deletion
            num_reads_other
  ).join("\t")

  # For each position in each reference
  Bio::FlatFile.foreach(options[:reference_fasta_file]) do |ref|
    ref_first_name = ref.definition.split(/\s/)[0]
    
    bioseq = ref.to_biosequence
    progress = ProgressBar.new(ref_first_name, bioseq.seq.length-3)
    
    bioseq.window_search(4,1) do |seq, position|
      progress.inc
      #next unless 54368 < position and 54380 > position
      
      
      log.debug "Seq: #{seq}, position #{position}" if log.debug?
      
      # If this position match a gag error,
      if options[:contexts_to_interrogate].include?(seq)
        log.debug "Apparently a context" if log.debug?
        region = "'#{ref_first_name}:#{position}-#{position+3}'"
        log.debug "Region: #{region}" if log.debug?

        pileup = sam.mpileup(:region => region){}
        
        if pileup.length==4 #happens when no sequences match the region being interrogated
          piles = Bio::DB::PileupIterator.new(pileup.join('')) #bit hacky, but eh.
          
          # For all reads that 
          pileups = piles.to_a
          reads = pileups.collect do |pile|
            pile.reads
          end.flatten.uniq
          
          # Look at the pileup for this region. How many are deletion, no deletion, or otherwise?
          # Output the contig name, position, num_with_deletion, num_without_deletion, num_otherwise
          num_reads_with_deletion = 0
          num_reads_without_deletion = 0
          num_reads_other = 0
          
          sequence_with_deletion = "#{seq[0..1]}#{seq[3]}"
          log.debug "deleted read to match to #{sequence_with_deletion}" if log.debug?
          reads.each do |read|
            current_read_sequence = read.sequence.gsub('*','')
            log.debug "read sequence: #{current_read_sequence}" if log.debug?
            
            if current_read_sequence == sequence_with_deletion              
              num_reads_with_deletion += 1
            elsif current_read_sequence == seq
              num_reads_without_deletion += 1
            else
              num_reads_other += 1
            end
          end
          
          if pileups[2].nil?
            $stderr.puts region
            $stderr.puts
            $stderr.puts pileup
            $stderr.puts pileups.inspect
          end
            
          
          puts [
            ref_first_name,
            position,
            seq,
            pileups[2].coverage,
            num_reads_with_deletion,
            num_reads_without_deletion,
            num_reads_other,
          ].join("\t")
        end
      end
    end
  end
  progress.finish
  
  
  
end #end if running as a script