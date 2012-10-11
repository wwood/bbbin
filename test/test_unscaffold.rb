require 'rspec'
require 'pp'
require 'open3'

# To run this test:
# $ rspec /path/to/test_script_being_tested.rb

# Assumes that the name of the file being tested is ../something.rb, and the name of this script is test_something.rb
$:.unshift File.join(File.dirname(__FILE__),'..')
script_under_test = File.basename(__FILE__).gsub(/^test_/,'')
require script_under_test
def assert_equal(e,o); o.should eq(e); end
path_to_script = File.join('..',script_under_test)


describe script_under_test do
  it 'should obvious test' do
    seqs = %w(>scaffold1 AANNTGTG)
    
    Open3.popen3(script_under_test) do |stdin, stdout, stderr|
      stdin.puts seqs.join("\n")
      stdin.close
      
      err = stderr.readlines
      raise err unless err == []
      answer = %w(>scaffold1_1of2_1_2 AA >scaffold1_2of2_5_8 TGTG).join("\n")+"\n"
      stdout.read.should eq(answer)
    end
  end
  
  it 'should handle multiple sequences' do
    seqs = %w(>scaffold1 AANNTGTG >s2 AAAAAAAA)
    
    Open3.popen3(script_under_test) do |stdin, stdout, stderr|
      stdin.puts seqs.join("\n")
      stdin.close
      
      err = stderr.readlines
      raise err unless err == []
      answer = %w(>scaffold1_1of2_1_2 AA >scaffold1_2of2_5_8 TGTG >s2_1of1_1_8 AAAAAAAA).join("\n")+"\n"
      stdout.read.should eq(answer)
    end
  end
end