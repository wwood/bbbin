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
  it 'should open3 test' do
    seqs =<<EOF
Methanoflorens_stordalmirensis_v4.3_scaffold3 822291  T 24  ,,,,.,...,,,..,...,,,,.,  QQQ\QQU\\\QQ\\RQQ\QQRB\Q
Methanoflorens_stordalmirensis_v4.3_scaffold3 822292  T 24  ,,,,.,...,,,..,...,,,,.,  QQQ^RQ\^^VQR^^RQQ^QRRB^Q
Methanoflorens_stordalmirensis_v4.3_scaffold3 822293  A 24  ,,,,.,...,,,..,...,,,,.,  &&&\&&\\\\&&\\'&&\&&'\\&
Methanoflorens_stordalmirensis_v4.3_scaffold3 822294  A 24  ggg,Gg...,gg..gGG.ggg,.g  !!!Z!!ZZZZ!!ZZ!!!Z!!!ZZ!
Methanoflorens_stordalmirensis_v4.3_scaffold3 822295  A 24  ,,,,.,...,,,..,...,,,,.,  !!!Z!!ZZZV!!ZZ!!!Z!!$RZ!
Methanoflorens_stordalmirensis_v4.3_scaffold3 822296  A 24  ,,,,.,...,,,..,...,,,,.,  !!!Z!!ZZZV!!ZZ!!!Z!!$RZ!
Methanoflorens_stordalmirensis_v4.3_scaffold3 822297  C 25  ,,,,.,...,,,..,...,,,,.,^]. TTThUTbcc^UUgiUTTiUUU^gUE
Methanoflorens_stordalmirensis_v4.3_scaffold3 822298  A 26  ,,,,.,...,,,..,...,,,,.,.^>,  TTThUTaccVUUdiUTTiUUU[eUUE
Methanoflorens_stordalmirensis_v4.3_scaffold3 822299  T 28  ,,,,.,...,,,..,...,,,,.,.,^>,^>,  __ghcdc_cHcbgh^ihic^cXecbUEE
Methanoflorens_stordalmirensis_v4.3_scaffold3 822300  G 28  ,,,,.,...,,,..,...,,,,.,.,,,  UUgicf[bbUcadhahhic^cWdce_UU
EOF
    seqs.gsub!(/[ ]+/, "\t")

    Open3.popen3(path_to_script) do |stdin, stdout, stderr|
      stdin.puts seqs
      stdin.close

      err = stderr.readlines
      err.should eq([]), err.join("")
      answer = 'TTAGAACATG'+"\n"
      stdout.read.should eq(answer)
    end
  end
end
