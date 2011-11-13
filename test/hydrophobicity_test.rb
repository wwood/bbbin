require "test/unit"
require 'hydrophobicity'

$:.unshift File.join(File.dirname(__FILE__),'..','ruby','lib')

# Doesn't work on ruby 1.8.6 - rails needs to be installed so the each_char method
# works
class HydrophobicityTest < Test::Unit::TestCase
  def test_hydrophobicity
    assert_equal [-0.4, -3.2, nil], Hydrophobicity.new.hydrophobicity_profile('GHX')
  end
  
  def test_acid_base
    assert_equal [nil, 10.54, nil], Hydrophobicity.new.profile('AKB', Hydrophobicity::PKA)
  end
end