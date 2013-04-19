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
  it 'when whole length is exactly a multiple of chunk size' do
    seqs = %w(>scaffold1 AANNTGTG)

    Open3.popen3("#{path_to_script} -c 4 -s 4") do |stdin, stdout, stderr|
      stdin.puts seqs.join("\n")
      stdin.close

      err = stderr.readlines
      err.should eq([]), err.join("")
      answer = %w(>scaffold1_1of2_1_4 AANN >scaffold1_2of2_5_8 TGTG).join("\n")+"\n"
      stdout.read.should eq(answer)
    end
  end

  it 'when whole length is not exactly a multiple of chunk size' do
    seqs = %w(>scaffold1 AANNTGTG)

    Open3.popen3("#{path_to_script} -c 4 -s 3") do |stdin, stdout, stderr|
      stdin.puts seqs.join("\n")
      stdin.close

      err = stderr.readlines
      err.should eq([]), err.join("")
      answer = %w(>scaffold1_1of3_1_4 AANN >scaffold1_2of3_4_7 NTGT >scaffold1_3of3_5_8 TGTG).join("\n")+"\n"
      stdout.read.should eq(answer)
    end
  end

  it 'bugs out, not quite sure how to solve this one' do
    seqs = %w(>scaffold1 AANNXTGTGXATGCAG)

    Open3.popen3("#{path_to_script} -c 10 -s 5") do |stdin, stdout, stderr|
      stdin.puts seqs.join("\n")
      stdin.close

      err = stderr.readlines
      err.should be_nil
      answer = %w(>scaffold1_1of3_1_10 AANNXTGTGX >scaffold1_2of3_6_15 TGTGXATGCA >scaffold1_3of3_7_16 GTGXATGCAG).join("\n")+"\n"
      stdout.read.should eq(answer)
    end
  end

  it 'should be ok printing the ends only' do
    seqs = %w(>scaffold1 AANNXTGTGXATGCAG >scaffold2 TAANNXTGTGXATGCAGT)

    Open3.popen3("#{path_to_script} --ends-only -c 3") do |stdin, stdout, stderr|
      stdin.puts seqs.join("\n")
      stdin.close

      err = stderr.readlines
      err.should eq([]), err
      answer = %w(>scaffold1_start AAN >scaffold1_end CAG
      >scaffold2_start TAA >scaffold2_end AGT).join("\n")+"\n"
      stdout.read.should eq(answer)
    end
  end
end
