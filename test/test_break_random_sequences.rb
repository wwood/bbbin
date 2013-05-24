require 'rspec'
require 'pp'
require 'systemu'
require 'tempfile'

# To run this test:
# $ rspec /path/to/test_script_being_tested.rb

# Assumes that the name of the file being tested is ../something.rb relative to the directory containing this test scripts, and the name of this tes script is test_something.rb
$:.unshift File.join(File.dirname(__FILE__),'..')
script_under_test = File.basename(__FILE__).gsub(/^test_/,'')
path_to_script = File.join(File.dirname(__FILE__),'..',script_under_test)


describe script_under_test do
  it 'should scripting test ok' do
    seqs = %w(>scaffold1_2of3_6_15 TGTGXATGCA >scaffold1_3of3_7_16 GTGXATGCAGAAAAAAAAAAAAAAAAAAAAAAAAAAA)

    Tempfile.open('break') do |tempfile|
      tempfile.puts seqs.join "\n"
      tempfile.close

      status, stdout, stderr = systemu "#{path_to_script} -n1 --min-length 6 #{tempfile.path}"
      raise stderr unless stderr == ""
      status.exitstatus.should eq(0)
      splits = stdout.split("\n")
      splits.length.should eq(6)
      splits[0].should eq('>scaffold1_2of3_6_15')
      splits[1].should eq('TGTGXATGCA')
      exp = '>scaffold1_3of3_7_16_break1_'
      splits[2][0...exp.length].should eq(exp)
      matches = splits[2].match(/>scaffold1_3of3_7_16_break1_(\d+)$/)

      raise unless matches
      s = matches[1].to_i
      s.should >= (6)
      s.should <= "GTGXATGCAGAAAAAAAAAAAAAAAAAAAAAAAAAAA".length-6
      exp = '>scaffold1_3of3_7_16_break2_'
      splits[4][0...exp.length].should eq(exp)
      "#{splits[3]}#{splits[5]}".should eq('GTGXATGCAGAAAAAAAAAAAAAAAAAAAAAAAAAAA')
    end
  end
end
