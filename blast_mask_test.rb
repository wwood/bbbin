require "test/unit"
require 'blast_mask'

class BlastMaskTest < Test::Unit::TestCase
  def test_mask_sequence1
    b = BlastHitArray.new
    assert_equal 'ABC', b.masked_sequence('ABC')
    b.push Hit.new(2,2)
    assert_equal 'AXC', b.masked_sequence('ABC')
    b.push Hit.new(1,3)
    assert_equal 'XXX', b.masked_sequence('ABC')
  end
end