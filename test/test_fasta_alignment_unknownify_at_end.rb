require 'rspec'
require 'pp'
require 'systemu'

# To run this test:
# $ rspec /path/to/test_script_being_tested.rb

# Assumes that the name of the file being tested is ../something.rb relative to the directory containing this test scripts, and the name of this tes script is test_something.rb
$:.unshift File.join(File.dirname(__FILE__),'..')
script_under_test = File.basename(__FILE__).gsub(/^test_/,'')
path_to_script = File.join(File.dirname(__FILE__),'..',script_under_test)


describe script_under_test do

  it 'should scripting test ok' do
    seqs = %w(>scaffold1 AAA).join "\n"

    status, stdout, stderr = systemu "#{path_to_script}", 'stdin' => seqs
    stderr.should eq("")
    status.exitstatus.should eq(0)
    stdout.chomp.should eq(seqs)
  end

  it 'should scripting test both sides' do
    seqs = %w(>scaffold1 --AAA-).join "\n"

    status, stdout, stderr = systemu "#{path_to_script}", 'stdin' => seqs
    stderr.should eq("")
    status.exitstatus.should eq(0)
    stdout.chomp.should eq(%w(>scaffold1 ..AAA.).join "\n")
  end
  it 'should scripting test LHS' do
    seqs = %w(>scaffold1 -----AAA).join "\n"

    status, stdout, stderr = systemu "#{path_to_script}", 'stdin' => seqs
    stderr.should eq("")
    status.exitstatus.should eq(0)
    stdout.chomp.should eq(%w(>scaffold1 .....AAA).join "\n")
  end
  it 'should scripting test RHS' do
    seqs = %w(>scaffold1 AGA--).join "\n"

    status, stdout, stderr = systemu "#{path_to_script}", 'stdin' => seqs
    stderr.should eq("")
    status.exitstatus.should eq(0)
    stdout.chomp.should eq(%w(>scaffold1 AGA..).join "\n")
  end
  it 'do deflines right' do
    seqs = ['>sca ffold me','--AA---'].join "\n"

    status, stdout, stderr = systemu "#{path_to_script}", 'stdin' => seqs
    stderr.should eq("")
    status.exitstatus.should eq(0)
    stdout.chomp.should eq(['>sca ffold me','..AA...'].join "\n")
  end
end
