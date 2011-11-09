require 'tempfile'
require 'rubygems'
require 'bio'

# Given 2 Bio:seq objects, blast them against each other and return the result.
module Bio
  class Blast

    class Bl2seq
      class Runner
        # Run a Bio::Seq (or Bio::FastaFormat) object against another. Assumes bl2seq is working correctly
        def bl2seq(seq1, seq2, options = {})
          optional_arguments = ''
          
          unless options.nil?
            if e = options[:evalue]
              optional_arguments += "-e #{e}"
            end
          end
          
          Tempfile.open('rubybl2seq') { |t1|  
            if seq1.kind_of?(Bio::FastaFormat)
              t1.puts seq1.to_s
            else
              t1.puts seq1.output(:fasta)
            end
            t1.close
            
            Tempfile.open('rubybl2seq') { |t2|
              if seq2.kind_of?(Bio::FastaFormat)
                t2.puts seq2.to_s
              else
                t2.puts seq2.output(:fasta)
              end
              
              t2.close
              
              Tempfile.open('rubybl2seqout') { |t3|  

                # Run the bl2seq. Assume protein blast for the moment
                cmd = "blastp #{optional_arguments} -query #{t1.path} -subject #{t2.path} -outfmt 6 -out #{t3.path}"
                ret = system cmd
                
              
                if !ret #Something went wrong
                  raise Exception, "Failed to run bl2seq: #{$?}"
                else
                  # Create the report from the output file
                  str = File.open(t3.path).read
                  
                  return Bio::Blast::Report.new(str)
                end
              }
            }
          }
        end
      end
    end
  end
end