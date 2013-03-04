require 'systemu'
require 'pp'
require 'open3'

# To run this test:
# $ rspec /path/to/test_script_being_tested.rb

# Assumes that the name of the file being tested is ../something.rb, and the name of this script is test_something.rb
$:.unshift File.join(File.dirname(__FILE__),'..')
script_under_test = File.basename(__FILE__).gsub(/^test_/,'')
#require script_under_test
def assert_equal(e,o); o.should eq(e); end
path_to_script = File.join(File.dirname(__FILE__),'..',script_under_test)

TEST_DATA_DIR = File.join(File.dirname(__FILE__),'..','testFiles','divvy_spectra')

describe script_under_test do
  it 'should do 1 protein hit' do
    test_file = "#{path_to_script} #{TEST_DATA_DIR}/single_protein.csv --trace error"
    status, stdout, stderr = systemu test_file

    stderr.should eq("")
    answer = "ID\tUnique spectra\tEstimated total spectra\tNormalised spectral count\tDescription\n"+
    ['Mstor_v4.3.2:1344','188','188.0','1.0','Methanoflorens_stordalmirensis_v4.3.2_01361 Methyl-coenzyme M reductase I subunit gamma '].join("\t")+"\n"
    stdout.should eq(answer)
  end

  it 'should do peptides that hit 2 proteins' do
    test_file = "#{path_to_script} #{TEST_DATA_DIR}/three_proteins.csv --trace error"
    status, stdout, stderr = systemu test_file

    stderr.should eq("")
    answer = "ID\tUnique spectra\tEstimated total spectra\tNormalised spectral count\tDescription\n"+
    ['Mstor_v4.3.2:1344','37','89.7075471698113','0.29509061569016876','Methanoflorens_stordalmirensis_v4.3.2_01361 Methyl-coenzyme M reductase I subunit gamma '+"\n"].join("\t")+
    ['eDeep20120820:eD1_8237_2','69','167.29245283018867','0.5503041211519364','methyl-coenzyme M reductase gamma subunit # pI:8.94 MW:27683 '+"\n"].join("\t")+
    ['eDeep20120820:eD1_1639_1','47','47.0','0.15460526315789475','chaperonin GroEL # pI:9.22 MW:10181 '+"\n"].join("\t")
    stdout.should eq(answer)
  end

  it 'should do peptides that hit more than 2 proteins' do
    test_file = "#{path_to_script} #{TEST_DATA_DIR}/multiply_mapped_spectra.csv --trace error"
    status, stdout, stderr = systemu test_file

    stderr.should eq("")
    answer = "ID\tUnique spectra\tEstimated total spectra\tNormalised spectral count\tDescription\n"+
    ['Mstor_v4.3.2:1344','37','89.7075471698113','0.24781090378400913','Methanoflorens_stordalmirensis_v4.3.2_01361 Methyl-coenzyme M reductase I subunit gamma '+"\n"].join("\t")+
    ['eDeep20120820:eD1_8237_2','69','167.29245283018867','0.46213384759720627','methyl-coenzyme M reductase gamma subunit # pI:8.94 MW:27683 '+"\n"].join("\t")+
    ['eDeep20120820:eD1_1639_1','47','47.0','0.1298342541436464','chaperonin GroEL # pI:9.22 MW:10181 '+"\n"].join("\t")+
    ['eDeep20120820:eD1_13975_5','9','9.157894736842106','0.02529805175923234','chaperonin GroEL # pI:6.58 MW:12040 '+"\n"].join("\t")+
    ['eDeep20120820:eD1_3396_1','38','38.666666666666664','0.10681399631675874','TGroEL # pI:9.70 MW:6451 '+"\n"].join("\t")+
    ['eDeep20120820:eD1_1494_8','10	10.175438596491228','0.028108946399147038','chaperonin GroEL # pI:4.93 MW:54101 '+"\n"].join("\t")

    stdout.should eq(answer)
  end

  it 'should ignore contaminants' do
    test_file = "#{path_to_script} #{TEST_DATA_DIR}/three_proteins_with_contaminant.csv --trace error"
    status, stdout, stderr = systemu test_file

    stderr.should eq("")
    answer = "ID\tUnique spectra\tEstimated total spectra\tNormalised spectral count\tDescription\n"+
    ['Mstor_v4.3.2:1344','37','89.7075471698113','0.29509061569016876','Methanoflorens_stordalmirensis_v4.3.2_01361 Methyl-coenzyme M reductase I subunit gamma '+"\n"].join("\t")+
    ['eDeep20120820:eD1_8237_2','69','167.29245283018867','0.5503041211519364','methyl-coenzyme M reductase gamma subunit # pI:8.94 MW:27683 '+"\n"].join("\t")+
    ['eDeep20120820:eD1_1639_1','47','47.0','0.15460526315789475','chaperonin GroEL # pI:9.22 MW:10181 '+"\n"].join("\t")
    stdout.should eq(answer)
  end

  it 'should do ok merging spectra' do
    test_file = "#{path_to_script} #{TEST_DATA_DIR}/three_proteins_meant_for_merge.csv --merge-proteins #{File.join(TEST_DATA_DIR,'merge_definition.csv')} --trace error"
    status, stdout, stderr = systemu test_file

    stderr.should eq("")
    answer = "ID\tUnique spectra\tEstimated total spectra\tNormalised spectral count\tDescription\n"+
    ['Mstor_v4.3.2:1344','37','89.7075471698113','0.29509061569016876','Methanoflorens_stordalmirensis_v4.3.2_01361 Methyl-coenzyme M reductase I subunit gamma '+"\n"].join("\t")+
    ['eDeep20120820:eD1_8237_2','69','167.29245283018867','0.5503041211519364','methyl-coenzyme M reductase gamma subunit # pI:8.94 MW:27683 '+"\n"].join("\t")+
    ['eDeep20120820:eD1_1639_1','47','47.0','0.15460526315789475','chaperonin GroEL # pI:9.22 MW:10181 '+"\n"].join("\t")
    stdout.should eq(answer)
  end
end
