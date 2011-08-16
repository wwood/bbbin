require 'test/unit'
require 'silica_sliding_window'

class SilicaSlidingWindowTest < Test::Unit::TestCase
  def test_window_is_hit?
    screener = Bio::SilaffinScreener.new
    assert_equal true, screener.window_is_hit?(Bio::Sequence::AA.new('SSSKK'))
    assert_equal true, screener.window_is_hit?(Bio::Sequence::AA.new('SSKKK'))
    assert_equal false, screener.window_is_hit?(Bio::Sequence::AA.new('KKKKK'))
    assert_equal false, screener.window_is_hit?(Bio::Sequence::AA.new('SSSSS'))
    assert_equal false, screener.window_is_hit?(Bio::Sequence::AA.new('SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS'))
    assert_equal true, screener.window_is_hit?(Bio::Sequence::AA.new('KKKKKKKKKK SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS'))
  end
  
  def test_hit?
    screener = Bio::SilaffinScreener.new
    
    assert_equal true, 
    screener.hit?(Bio::Sequence::AA.new('KKKKKKKKKK SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS KKKKKKK')),
    '100 amino acids 10% K, 90% S'
    
    assert_equal true, 
    screener.hit?(Bio::Sequence::AA.new('SKKKKKKKKK SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS K KKKK KKKKKKKKKK KKKKKKKKKK KKKKKKKKKK')),
    '101 amino acids 10% K, 90% S 100 stretch'
    
    assert_equal true, 
    screener.hit?(Bio::Sequence::AA.new('SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS KKKKKKKKKS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS KKKKKKKKKK SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS KKKKKKKKKS')),
    '120 amino acids 10% K, 90% S in a longer than 100aa stretch'
    
    assert_equal false, 
    screener.hit?(Bio::Sequence::AA.new('KKKKKKKKKK KKKKKKKKKK KKKKKKKKKK KKKKKKKKKK KKKKKKKKKK KKKKKKKKKK KKKKKKKKKK SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS')),
    'no signal peptide'
    
    assert_equal false, 
    screener.hit?(Bio::Sequence::AA.new('KKKKKKKKKK SSSSSSSSSS SSSSSSSSSS SSSSSAAAAA AAAAAAAAAA SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS SSSSSSSSSS')),
    'when signal peptide is removed, there is no longer a satisfying window'
  end
end