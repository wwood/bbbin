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
path_to_script = File.join(File.dirname(__FILE__),'..',script_under_test)



describe script_under_test do
  
  it 'should test forward first' do
    Open3.popen3(path_to_script) do |stdin, stdout, stderr|
      input = '@E1D3_10_1/1 reference=NC_008278 position=4711992..4712096 errors=76- description="Frankia alni ACN14a, complete genome."
CCCTCGCCGGTCACCGCCCGCCGCAGGAAGCCGCCCACGCCCTGGCCCTTGTAAGCGAAGTCGAGCCGCCCCCGGAGGCCACCATCGACCCCTGCTTGGCGTAG
+
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
      stdin.puts input
      stdin.close
      
      err = stderr.readlines
      err.should eq([]), err.join("")
      answer = input+"\n"
      stdout.read.should eq(answer)
    end
  end

  it 'should test forward second' do
    Open3.popen3(path_to_script) do |stdin, stdout, stderr|
      input = '@E1D3_10_1/2 reference=NC_008278 position=complement(4712355..4712459) errors=95+G description="Frankia alni ACN14a, complete genome."
TCGCGTGCGGTGCCATCGCGTGCCGTGCCATCGCGTGCGGTGCCGTCGACGCGGAGCCGGGGACGCCGTCCTGCCGGCTTGCCGGGCCGAGGCCTGCCTTCCGACG
+
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXSXXXXXXXXXX'
      stdin.puts input
      stdin.close
      
      err = stderr.readlines
      err.should eq([]), err.join("")
      answer = input+"\n"
      stdout.read.should eq(answer)
    end
  end
  
  
  
  it 'should test reverse first' do
    Open3.popen3(path_to_script) do |stdin, stdout, stderr|
      input = '@E1D3_10_2/1 reference=NC_013715 position=complement(27203..27307) errors=59%T description="Rothia mucilaginosa DY-18, complete genome."
TGATTTCGTCGCGGCGGGTAATGGCGGAGCACATTTCCTTCACGCGGGCATAGGTCTCTGGGTCGGTGGTTGCCTTCGCTTCAGCCTGCAGCGCCTCGAGGGCTG
+
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXSXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
      stdin.puts input
      stdin.close
      
      err = stderr.readlines
      err.should eq([]), err.join("")
      
      expected = '@E1D3_10_2/1 reference=NC_013715 position=complement(27203..27307) errors=59%T description="Rothia mucilaginosa DY-18, complete genome."
CAGCCCTCGAGGCGCTGCAGGCTGAAGCGAAGGCAACCACCGACCCAGAGACCTATGCCCGCGTGAAGGAAATGTGCTCCGCCATTACCCGCCGCGACGAAATCA
+
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXSXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
      answer = expected+"\n"
      stdout.read.should eq(answer)
    end
  end
  
  
  it 'should test reverse second' do
    Open3.popen3(path_to_script) do |stdin, stdout, stderr|
      input = '@E1D3_10_2/2 reference=NC_013715 position=27604..27708 errors=105%T description="Rothia mucilaginosa DY-18, complete genome."
GCAATCCCCGTGGTGGAAAACGCGCCCGCCTCGCCCCTGGTCACTAACGAGGTAACCCCCAAGTCCTAGCCGGTACGCCCGTGCGCGAACACTACGTTTCTGCCT
+
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXS'
      stdin.puts input
      stdin.close
      
      err = stderr.readlines
      err.should eq([]), err.join("")
      
      expected = '@E1D3_10_2/2 reference=NC_013715 position=27604..27708 errors=105%T description="Rothia mucilaginosa DY-18, complete genome."
AGGCAGAAACGTAGTGTTCGCGCACGGGCGTACCGGCTAGGACTTGGGGGTTACCTCGTTAGTGACCAGGGGCGAGGCGGGCGCGTTTTCCACCACGGGGATTGC
+
SXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
      answer = expected+"\n"
      stdout.read.should eq(answer)
    end
  end
  
  
  

end