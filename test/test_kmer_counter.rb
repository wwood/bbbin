#!/usr/bin/env ruby

$:.unshift File.join(File.dirname(__FILE__),'..')

require 'test/unit'
require 'kmer_counter'
require 'tempfile'

class KmerCounterTest < Test::Unit::TestCase
  def test_lowest_lexigraphical_form
    assert_equal Bio::Sequence::NA.new('AA'), Bio::Sequence::NA.new('AA').lowest_lexigraphical_form
    assert_equal Bio::Sequence::NA.new('AA'), Bio::Sequence::NA.new('TT').lowest_lexigraphical_form
    assert_equal Bio::Sequence::NA.new('AG'), Bio::Sequence::NA.new('CT').lowest_lexigraphical_form
  end
  
  def test_empty_kmer_hash
    answer = {}; %w(A C).each{|k| answer[k] = 0}
    assert_equal answer, Bio::Sequence::Kmer.empty_kmer_hash(1)
    assert_equal 136, Bio::Sequence::Kmer.empty_kmer_hash.length
  end
  
  def script_path
    File.join(File.dirname(__FILE__),'..','kmer_counter.rb')
  end
  
  def test_running1
    Tempfile.open('one') do |tempfile|
      tempfile.puts '>one'
      tempfile.puts 'ACAGT'
      tempfile.close

      assert_equal "ID\tA\tC\none_0\t0.6\t0.4\n", `#{script_path} -w 5 -k 1 #{tempfile.path}`
    end
  end
end