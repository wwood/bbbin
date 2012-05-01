
#$:.unshift File.join(File.dirname(__FILE__),'..')

require 'test/unit'
require 'tmpdir'

class BlastBySplitsTest < Test::Unit::TestCase
  def setup
    @blast_by_splits = "#{File.expand_path File.join(File.dirname(__FILE__),'..','blast_by_splits.rb')}"
  end
  
  def test_evalue
    acp = ">acp\nMKILLLCIIFLYYVNAFKNTQKDGVSLQILKKKRSNQVNFLNRKNDYNLIKNKNPSSSLKSTFDDIKKIISKQLSVEEDK"
  
    # no idea if this stuff works on platforms other than linux
    Dir.mktmpdir do |dir|
      Dir.chdir dir
      f = File.open('acp.fa', 'w')
      f.puts acp
      f.close
      
      raise unless `makeblastdb -in acp.fa -dbtype prot`
      raise unless `#{@blast_by_splits} -a 2 -o out -i acp.fa -d acp.fa`
      assert_equal "acp\tacp\t100.00\t80\t0\t0\t1\t80\t1\t80\t1e-43\t 154\n", File.open('out').read

      raise unless `#{@blast_by_splits} -a 2 -o out -i acp.fa -d acp.fa -e 1e-50`
      assert_equal "", File.open('out').read
    end
  end
end
