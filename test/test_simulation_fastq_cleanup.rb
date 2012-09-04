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
    seqs = %w(>scaffold1 AANNTGTG)
    Tempfile.open('fastq') do |f|
      f.puts <<EOF
@1/1 reference=NC_003454 position=complement(646364..646468) errors=80%A description="Fusobacterium nucleatum subsp. nucleatum ATCC 25586, complete genome."
TATCTTCCACTAAAATTTCTTCTACACAATCTTGAATTAAACTAATATTTTCAGTTTTTTCAAGTTTTTCTCTCATTTTACTTCTATATTTATATTTATCTGCCT
+
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXSXXXXXXXXXXXXXXXXXXXXXXXXX
@1/2 reference=NC_003454 position=646804..646908 errors=94%A description="Fusobacterium nucleatum subsp. nucleatum ATCC 25586, complete genome."
AAGAATAATACTGTTCCTACTTGGCTTACATATACTTCAGACAAAACTATTGAAGTTATAAAAGAAATGATGAAATTTTCTCCAATAGTTAGTAGAATGGTAAAT
+
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXSXXXXXXXXXXX
@2/1 reference=NC_015160 position=3262693..3262797 errors=98%G description="Odoribacter splanchnicus DSM 20712 chromosome, complete genome."
CCCTCCTACAAAGGTCAGATTTTGGTTCCTACCTATCCTATGATCGGTAATT
+
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
@2/2 reference=NC_015160 position=3262693..3262797 errors=98%G description="Odoribacter splanchnicus DSM 20712 chromosome, complete genome."
CCCTCCTACAAAGGTCAGATTTTGGTTCCTACCTATCCTATGATCGGTAATT
+
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
EOF
      f.close

      Open3.popen3("#{path_to_script} -n E2F #{f.path} #{f.path}") do |stdin, stdout, stderr|
        err = stderr.readlines
        err.should eq([]), err.join("")
        answer =<<EOF
@1/1 E2F
TATCTTCCACTAAAATTTCTTCTACACAATCTTGAATTAAACTAATATTTTCAGTTTTTTCAAGTTTTTCTCTCATTTTACTTCTATATTTATATTTATC
+
fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffTffffffffffffffffffff
@1/2 E2F
AAGAATAATACTGTTCCTACTTGGCTTACATATACTTCAGACAAAACTATTGAAGTTATAAAAGAAATGATGAAATTTTCTCCAATAGTTAGTAGAATGG
+
fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffTffffff
@2/1 E2F
CCCTCCTACAAAGGTCAGATTTTGGTTCCTACCTATCCTATGATCGGTAATT
+
ffffffffffffffffffffffffffffffffffffffffffffffffffff
@2/2 E2F
CCCTCCTACAAAGGTCAGATTTTGGTTCCTACCTATCCTATGATCGGTAATT
+
ffffffffffffffffffffffffffffffffffffffffffffffffffff
@3/1 E2F
TATCTTCCACTAAAATTTCTTCTACACAATCTTGAATTAAACTAATATTTTCAGTTTTTTCAAGTTTTTCTCTCATTTTACTTCTATATTTATATTTATC
+
fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffTffffffffffffffffffff
@3/2 E2F
AAGAATAATACTGTTCCTACTTGGCTTACATATACTTCAGACAAAACTATTGAAGTTATAAAAGAAATGATGAAATTTTCTCCAATAGTTAGTAGAATGG
+
fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffTffffff
@4/1 E2F
CCCTCCTACAAAGGTCAGATTTTGGTTCCTACCTATCCTATGATCGGTAATT
+
ffffffffffffffffffffffffffffffffffffffffffffffffffff
@4/2 E2F
CCCTCCTACAAAGGTCAGATTTTGGTTCCTACCTATCCTATGATCGGTAATT
+
ffffffffffffffffffffffffffffffffffffffffffffffffffff
EOF

        stdout.read.should eq(answer)
      end
    end
  end
end