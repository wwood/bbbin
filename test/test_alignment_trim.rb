require 'rspec'

# Assumes that the name of the file being tested is ../something.rb, and the name of this script is test_something.rb
$:.unshift File.join(File.dirname(__FILE__),'..')
script_under_test = File.basename(__FILE__).gsub(/^test_/,'')
require script_under_test
def assert_equal(e,o); o.should eq(e); end

describe script_under_test do
  it 'should obvious test' do
    seqs = %w(a a a)
    trimmed = Bio::Alignment.new(seqs).trim_uninformative_columns
    trimmed[0].size.should eq(0) #"all seqs should be empty since they are all uninformative"
    assert_equal 3, trimmed.size #"should still have 3 seqs in it"
  end
  
  it 'retain informative ones' do
    seqs = %w(atg agt ttc)
    trimmed = Bio::Alignment.new(seqs).trim_uninformative_columns
    assert_equal 'g', trimmed[0]
    assert_equal 't', trimmed[1]
  end
  
  it 'remove empty columns' do
    seqs = %w(at--g ag--t tt--c)
    trimmed = Bio::Alignment.new(seqs).trim_uninformative_columns
    assert_equal 'g', trimmed[0]
  end
  
  it 'remove almost empty columns' do
    seqs = %w(at--g ag-at tt--c)
    trimmed = Bio::Alignment.new(seqs).trim_uninformative_columns
    assert_equal 'g', trimmed[0]
  end
  
  it 'retain more than one column' do
    seqs = %w(atagy agttk ttgcr)
    trimmed = Bio::Alignment.new(seqs).trim_uninformative_columns
    assert_equal 'agy', trimmed[0]
  end
  
  it "chuck out columns with too many empty chars" do
    seqs = %w(atg g-t tac)
    trimmed = Bio::Alignment.new(seqs).trim_uninformative_columns(:maximum_missing_entries => 1.0)
    assert_equal 'atg', trimmed[0]
    trimmed = Bio::Alignment.new(seqs).trim_uninformative_columns(:maximum_missing_entries => 0.0)
    assert_equal 'ag', trimmed[0]
    trimmed = Bio::Alignment.new(seqs).trim_uninformative_columns(:maximum_missing_entries => 0.5)
    assert_equal 'atg', trimmed[0]
    trimmed = Bio::Alignment.new(seqs).trim_uninformative_columns(:maximum_missing_entries => 0.1)
    assert_equal 'ag', trimmed[0]
  end
  
  it "should work with non-standard missing characters" do
    seqs = %w(atg g-t tac)
    trimmed = Bio::Alignment.new(seqs).trim_uninformative_columns(:maximum_missing_entries => 0.0, :missing_character => '?')
    assert_equal 'atg', trimmed[0]
  end
  
  it 'test_running1' do
    Tempfile.open('one') do |tempfile|
      tempfile.puts '>1'
      tempfile.puts 'atg'
      tempfile.puts '>2'
      tempfile.puts 'g-t'
      tempfile.puts '>3'
      tempfile.puts 'tac'
      tempfile.close

      assert_equal %w(>1 atg >2 g-t >3 tac).join("\n")+"\n", `#{script_under_test} #{tempfile.path}`
      assert_equal %w(>1 ag >2 gt >3 tc).join("\n")+"\n", `#{script_under_test} -m 0.1 #{tempfile.path}`
    end
  end
  
end