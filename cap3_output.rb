#!/usr/bin/ruby

# A script for parsing the stdout of CAP3 assembly runs. It gives
# the links between the input traces, as well as the assembly details,
# just in string form
module Bio
  class Cap3
    class Output
      attr_reader :hierarchy, :assembly
  
  
      # parse an output file into a hash of traces and assemblies
      def parse_file(filename)
        parse(File.open(filename).read)
      end
      
      def parse(stream)
        lines = stream.split("\n")
        cur_contig = nil
        @hierarchy = {}
        i = 0;
    
    
        # Create the hashes
        lines.each do |line|
          i += 1
          #      line.strip!
      
          # if *** type line, change the contig name we are working with
          if match = line.match('^\*+ (.+) \*+$')
            cur_contig = match[1]
            next
          end
      
          # if no current contig or just plain nothn, then skip this line
          if !cur_contig or line.strip.empty?
            next
          end
      
          # No more links, go to next stage
          if line.match('DETAILED DISPLAY OF CONTIGS')
            break
          end
      
          # line is a link. should be either like these:
          # CHR4524-1+
          #              CHR1739+ is in CHR4524-1+
          if match = line.match('\s+(.+?)[\+\-] is in ')
            @hierarchy[cur_contig].push match[1]
            next
          elsif match = line.match('^(.+)[\-\+]$')
            if @hierarchy[cur_contig] #if no entry, need to insert an array, rather than pushing to an already created array
              @hierarchy[cur_contig].push match[1]
            else
              @hierarchy[cur_contig] = [match[1]]
            end
          else
            raise Exception, "Badly Parsed File, on this line: #{line}"
          end
        end
    
    
    
        # Record each of the assembly details
        @assembly = {}
        cur_contig = nil
    
        lines[i..lines.length-1].each do |line|
          line.chomp!      
          if match = line.match('^\*+ (.+) \*+$')
            cur_contig = match[1]
            @assembly[cur_contig] = []
          else
            @assembly[cur_contig].push line
          end
        end
        
        return self
      end
      
      def parse_concatenated(stream)
        lines = stream.split("\n")
        cur_contig = nil
        @hierarchy = {}
        @assembly = {}
        hierarchy_mode = true
        i = 0;
    
    
        # Create the hashes
        lines.each do |line|
          i += 1
      
          # if *** type line, change the contig name we are working with
          if match = line.match('^\*+ (.+) \*+$')
            cur_contig = match[1]
            if @hierarchy[cur_contig]
              hierarchy_mode = false
            else
              hierarchy_mode = true
            end
            next
          end

          # if no current contig or just plain nothn, then skip this line
          if !cur_contig or line.strip.empty?
            next
          end
          
          if hierarchy_mode
            # line is a link. should be either like these:
            # CHR4524-1+
            #              CHR1739+ is in CHR4524-1+
            if match = line.match('\s+(.+?)[\+\-] is in ')
              @hierarchy[cur_contig].push match[1]
              next
            elsif match = line.match('^(.+)[\-\+]$')
              if @hierarchy[cur_contig] #if no entry, need to insert an array, rather than pushing to an already created array
                @hierarchy[cur_contig].push match[1]
              else
                @hierarchy[cur_contig] = [match[1]]
              end
            elsif line.match(/DETAILED DISPLAY OF CONTIGS/)
              #ignore
            else
              raise Exception, "Badly Parsed File, on this line: #{line}"
            end
            
            
          else #assembly mode
            line.chomp!
            if  @assembly[cur_contig]
              @assembly[cur_contig].push line
            else
              @assembly[cur_contig] = [line]
            end
          end
        end
    
    
    
        # Record each of the assembly details
        return self
      end
    end
  end
end


# If called as a script, print out the contig names, the number of traces
# and then a comma separated list of the trace names
if $0 == __FILE__
  p = Bio::Cap3::Output.new
  p.parse_concatenated(ARGF.read)
  
  p.hierarchy.to_a.sort{|a,b| a[0]<=>b[0]}.each do |contig, traces|
    puts [
      contig,
      traces.length,
      traces.join(', ')
    ].join("\t")
  end
end
