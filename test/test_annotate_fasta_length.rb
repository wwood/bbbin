require 'rspec'
require 'pp'
require 'open3'
require 'tempfile'
# To run this test:
# $ rspec /path/to/test_script_being_tested.rb

# Assumes that the name of the file being tested is ../something.rb, and the name of this script is test_something.rb
$:.unshift File.join(File.dirname(__FILE__),'..')
script_under_test = File.basename(__FILE__).gsub(/^test_/,'')
require script_under_test
def assert_equal(e,o); o.should eq(e); end
path_to_script = File.join(File.dirname(__FILE__),'..',script_under_test)



describe script_under_test do
  it 'should open3 test' do
    seqs = [
      '>one seq',
      'ATG',
      '>two',
      'A'*80,
      'T'*40,
    ]
    
    Open3.popen3(path_to_script) do |stdin, stdout, stderr|
      stdin.puts seqs.join("\n")
      stdin.close
      
      err = stderr.readlines
      err.should eq([]), err.join("")
      answer = [
        '>one_len3 seq',
        'ATG',
        '>two_len120',
        'A'*80,
        'T'*40,
      ].join("\n")+"\n"
      stdout.read.should eq(answer)
    end
  end
  
  it 'should accept a fasta file as an argument' do
    seqs = [
      '>one seq',
      'ATG',
      '>two',
      'A'*80,
      'T'*40,
    ]
    
    Tempfile.open('test') do |tempfile|
      tempfile.puts seqs.join "\n"
      tempfile.close
      
      Open3.popen3(path_to_script+" #{tempfile.path}") do |stdin, stdout, stderr|
        stdin.close
        
        err = stderr.readlines
        err.should eq([]), err.join("")
        answer = [
          '>one_len3 seq',
          'ATG',
          '>two_len120',
          'A'*80,
          'T'*40,
        ].join("\n")+"\n"
        stdout.read.should eq(answer)
      end
      
    end

  end
end