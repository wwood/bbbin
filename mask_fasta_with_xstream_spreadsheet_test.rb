require "test/unit"
require 'rubygems'
require 'construct'

class MaskFastaWithXstreamSpreadsheetTest < Test::Unit::TestCase
  include Construct::Helpers
  
  def setup
    @test_data_dir = File.join('/home/ben/bin/testFiles') 
  end
  
  def test_three_with_long_ids
    # Assumes the script under test is in the path
    within_construct do |construct|
      construct.directory 'alice/rabbithole' do |dir|
        system "mask_fasta_with_xstream_spreadsheet.rb #{File.join(@test_data_dir,'mask_fasta_with_xstream_spreadsheet.fa')} #{File.join(@test_data_dir,'XSTREAM__i0.7_g3_m1_e3.0_chart.xls')} >out"
        assert_equal '>PF123 456 | blah
XXXXXXXXXXXXXXXXXXXXXXXXBDSSXXXXX
>PF123 dsa456 | blah
ASDEDFA
>PF123 456vfd | blah
AXXXXXXBDSSXXXXX
',
          File.read('out')
      end
    end
  end
  
end