require "test/unit"

$:.unshift File.join(File.dirname(__FILE__),'..')
require 'gly_sliding_window'

class SequenceWindowDescriptorTest < Test::Unit::TestCase
  def test_first
    desc = Bio::SequenceWindowDescriptor.new
    desc.calculate(Bio::Sequence::AA.new('GGY'), 2)
    assert_equal 2, desc.maximum_counts[:gly]
    assert_equal 'GG', desc.maximum_sequences[:gly]
  end

  def test_not_first
    desc = Bio::SequenceWindowDescriptor.new
    desc.calculate(Bio::Sequence::AA.new('GGYYYGGYGY'), 4)
    assert_equal 3, desc.maximum_counts[:gly]
    assert_equal 'GGYG', desc.maximum_sequences[:gly]
  end

  def test_not_big_enough
    desc = Bio::SequenceWindowDescriptor.new
    desc.calculate(Bio::Sequence::AA.new('GG'), 4)
    assert_equal 0, desc.maximum_counts[:gly]
    assert_equal '', desc.maximum_sequences[:gly]
  end
end