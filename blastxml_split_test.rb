# To change this template, choose Tools | Templates
# and open the template in the editor.

$:.unshift File.join(File.dirname(__FILE__),'..','lib')

require 'test/unit'

class BlastxmlSplitTest < Test::Unit::TestCase
  def test_multifasta
    # parse a multi-fasta blastxml output directly, and then parse each of the splits, and compare minus the headers
    assert File.exists?('1.xml') == false
    assert File.exists?('2.xml') == false
    assert File.exists?('blastxml_split_test_direct') == false

    `blastxml_to_tab.rb testFiles/multifasta_blast.xml |grep -v 'Query Definition' >blastxml_split_test_direct`
    `blastxml_split.rb testFiles/multifasta_blast.xml`
    assert File.exists?('1.xml')
    assert File.exists?('2.xml')
    assert File.exists?('3.xml') == false
    `blastxml_to_tab.rb 1.xml >1.csv`
    `blastxml_to_tab.rb 2.xml >2.csv`
    `cat 1.csv 2.csv |grep -v 'Query Definition' |diff - blastxml_split_test_direct >diff`
    assert_equal '', File.open('diff').read

    `rm 1.xml`
    `rm 2.xml`
    `rm 1.csv`
    `rm 2.csv`
    `rm diff`
    `rm blastxml_split_test_direct`
  end
end
