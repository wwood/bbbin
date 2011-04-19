#!/usr/bin/env/ruby

require 'test/unit'
require 'primer_design_extras'

class PrimerExtrasTest < Test::Unit::TestCase
  # def test_oligo_designer
  #   o = OligoDesigner.new
  #   assert_equal 63.627401, o.melting_temperature('TGTGGCAGCTGCTTAGTGTAGCGCG')
  #   assert_equal 'TGTGGCAGCTGCT', o.just_below('TGTGGCAGCTGCTTAGTGTAGCGCG', 44.5)
  #   assert_equal 'TGTGGCAGCTGCTTAGTGTAGCGCG', o.just_below('TGTGGCAGCTGCTTAGTGTAGCGCG', 150)
  #   assert_equal 'TGTGGCAGCTGCTTAGTGTAGCGC', o.just_below('TGTGGCAGCTGCTTAGTGTAGCGCG', 61.8)
  # end

  # def test_oligo_possible_oligos_ordered_by_temperature_difference
  #   # ben@uyen:~/bin$ oligotm -tp 1 -sc 1 TGTGGCAGCTGCTTAGTGTAGCGCG
  #   # 63.627401
  #   # ben@uyen:~/bin$ oligotm -tp 1 -sc 1 TGTGGCAGCTGCTTAGTGTAGCGC
  #   # 61.784774
  #   # ben@uyen:~/bin$ oligotm -tp 1 -sc 1 TGTGGCAGCTGCTTAGTGTAGCG
  #   # 59.514664
  #   # ben@uyen:~/bin$ oligotm -tp 1 -sc 1 TGTGGCAGCTGCTTAGTGTAGC
  #   # 57.262267
  #   # ben@uyen:~/bin$ oligotm -tp 1 -sc 1 TGTGGCAGCTGCTTAGTGTAG
  #   # 54.517682
  #   # ben@uyen:~/bin$ oligotm -tp 1 -sc 1 TGTGGCAGCTGCTTAGTGTA
  #   # 53.463947
  #   o = OligoDesigner.new
  #   # test no possible
    # assert_equal [], o.possible_oligos_ordered_by_temperature_difference('TGTGGCAGCTGCTTAGTGTAGCGCG',64,66,68,2)
    # # test one possible
    # assert_equal %w(TGTGGCAGCTGCTTAGTGTAGCGC), o.order('TGTGGCAGCTGCTTAGTGTAGCGCG',60,61,62,2)
    # # test three, with ordering
    # assert_equal %w(
    # TGTGGCAGCTGCTTAGTGTAGCGCG
    # TGTGGCAGCTGCTTAGTGTAGCGC
    # TGTGGCAGCTGCTTAGTGTAGCG
    # ), o.order('TGTGGCAGCTGCTTAGTGTAGCGCG',58,63.5,64,2), 'order 1'
    # # test three, with different ordering
    # assert_equal %w(
    # TGTGGCAGCTGCTTAGTGTAGCG
    # TGTGGCAGCTGCTTAGTGTAGCGC
    # TGTGGCAGCTGCTTAGTGTAGCGCG
    # ), o.order('TGTGGCAGCTGCTTAGTGTAGCGCG',58,58,64,2), 'order 2'
 #   # test two, with GC clamp
  #   assert_equal %w(
  #   TGTGGCAGCTGCTTAGTGTAGCGCG
  #   TGTGGCAGCTGCTTAGTGTAGCGC
  #   ), o.order('TGTGGCAGCTGCTTAGTGTAGCGCG',58,63.5,64,4), 'gc clamp 4'
  # end

  def test_boulder_io
    assert_equal "three=four\none=two\n=\n", BoulderIO::Record.new({'one'=>'two','three'=>'four'}).to_s
    assert BoulderIO::Record.new({'one'=>'two','three'=>'four'}).is_a?(BoulderIO::Record)
    assert 'four', BoulderIO::Record.new({'one'=>'two','three'=>'four'})['three']
    bio = BoulderIO::Record.new({'one'=>'two','three'=>'four'})
    bio['you say 911'] = 'I say 000' #emergency services in aussie-land is 000
    # it's much better I agree
    assert 'I say 000', bio['you say 911']
  end

  def test_boulder_io_reading
    # The test file:
    # PRIMER_SEQUENCE_ID=test1
    # SEQUENCE=GACTGATCGATGCTAGCTACGATCGATCGATGCATGCTAGCTAGCTAGCTGCTAGC
    #=
    count = 0
    BoulderIO.open('testFiles/boulderio_example.txt').each do |record|
      assert record.is_a?(BoulderIO::Record)
      assert_equal 'test1', record['PRIMER_SEQUENCE_ID']
      assert_equal 'GACTGATCGATGCTAGCTACGATCGATCGATGCATGCTAGCTAGCTAGCTGCTAGC', record['SEQUENCE']
      assert_equal 2, record.attributes.length
      count += 1
    end
    assert_equal 1, count
  end
  
  def test_primer3_output_result
  	res = Primer3Result.create_from_primer3_output_filename('testFiles/primer3output_negative.txt')
  	assert_equal false, res.yeh?
  	
  	res = Primer3Result.create_from_primer3_output_filename('testFiles/primer3output_positive.txt')
  	assert_equal true, res.yeh?
  	assert_equal '59.997', res['PRIMER_LEFT_TM']
  end
  
  def test_PossiblyRestrictedNucleotideSequence_unique_within_region?
  	e1 = EnzymeCut.new('NeckoR2D2',3,5)
  	pr = 	PossiblyRestrictedNucleotideSequence.new([e1])
  	assert_equal true, pr.unique_within_region?(e1,1,10), "same enzyme"
  	assert_equal false, pr.unique_within_region?(EnzymeCut.new('NeckoR2D2',4,6),1,10), "not unique"
  	assert_equal true, pr.unique_within_region?(EnzymeCut.new('NeckR2E2',3,4),1,10), "different enzyme"
  	assert_equal true, pr.unique_within_region?(EnzymeCut.new('NeckoR2D2',13,14),1,10), "outside bounds upper"
  	assert_equal true, pr.unique_within_region?(EnzymeCut.new('NeckoR2D2',1,2),5,10), "outside bounds lower"
  	assert_equal false, pr.unique_within_region?(EnzymeCut.new('NeckoR2D2',7,8),1,10), "another fail"
  end
end