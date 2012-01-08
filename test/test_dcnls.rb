#!/usr/bin/env ruby

$:.unshift File.join(File.dirname(__FILE__),'..')

require 'test/unit'
require 'dcnls'

class DCNLSTest < Test::Unit::TestCase
  def test_single
    assert_equal [2,3],
    Bio::DCNLS.new.predictions(Bio::Alignment.new(%w(
    GGKKKKGG
    )))
    
    assert_equal [],
    Bio::DCNLS.new.predictions(Bio::Alignment.new(%w(
    GGKGKGG
    )))
  end
  
  def test_msa
    assert_equal [],
    Bio::DCNLS.new.predictions(Bio::Alignment.new(%w(
    GGKKKKGG
    GGRGGRRG
    )))
    
    assert_equal [2,3],
    Bio::DCNLS.new.predictions(Bio::Alignment.new(%w(
    GGKKKKGG
    GGRRRRRG
    )))
  end
  
  def test_tweak_required_number_of_basic_residues
    assert_equal [],
    Bio::DCNLS.new.predictions(Bio::Alignment.new(%w(
    GGKKKKGG
    )),
    :required_number_of_basic_residues => 5
    )
  end
end