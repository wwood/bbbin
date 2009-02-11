$:.unshift File.join(File.dirname(__FILE__),'..','lib')

require 'test/unit'
require 'rubygems'
require 'bio'
require 'orf_finder'
include Orf

class OrfFinderTest < Test::Unit::TestCase
  def test_simples
    finder = OrfFinder.new
    
    # simple as
    threads = finder.generate_longest_orfs('ATGTAG') #start stop
    assert_equal 3, threads.length
    assert_kind_of OrfThread, threads[1]
    assert_equal 1, threads[0].length
    t = threads[0]
    assert_equal 0, t[0].start
    assert_equal 5, t[0].stop
    assert_equal 0, threads[1].length
    assert_equal 0, threads[2].length
    
    # off 1 frame
    threads = finder.generate_longest_orfs('AATGTAG') #start stop
    assert_equal 3, threads.length
    assert_kind_of OrfThread, threads[1]
    assert_equal 1, threads[1].length
    t = threads[1]
    assert_equal 1, t[0].start
    assert_equal 6, t[0].stop
    assert_equal 0, threads[0].length
    assert_equal 0, threads[2].length
    
    # has an partial frame at the end
    threads = finder.generate_longest_orfs('ATGAAATAGATGAAA')
    assert_equal 3, threads.length
    assert_kind_of OrfThread, threads[1]
    assert_equal 2, threads[0].length
    assert_equal 1, threads[1].length
    assert_equal 0, threads[2].length
    t = threads[0]
    assert_equal 0, t[0].start
    assert_equal 8, t[0].stop
    assert_equal false, t[0].fragment?, t[0].aa_sequence
    assert_equal 9, t[1].start
    assert_equal 14, t[1].stop
    assert_equal true, t[1].fragment?
    
    # has partial frame at start and end in the third frame
    #irb(main):004:0> Bio::Sequence::NA.new('AAAAATAGATGAAATAGATGAAT').translate(3)
    #=> "K*MK*MN"
    #irb(main):005:0> Bio::Sequence::NA.new('AAAAATAGATGAAATAGATGAAT').translate(1)
    #=> "KNR*NR*"
    #irb(main):006:0> Bio::Sequence::NA.new('AAAAATAGATGAAATAGATGAAT').translate(2)
    #=> "KIDEIDE"
    threads = finder.generate_longest_orfs('AAAAATAGATGAAATAGATGAAT') #translate(3) => "K*MK*MN"
    assert_equal 1, threads[0].length
    assert_equal 0, threads[1].length
    assert_equal 3, threads[2].length
    t = threads[0]
    assert_equal 0, t[0].start
    assert_equal 11, t[0].stop
    assert_equal true, t[0].fragment?
    t = threads[2]
    assert_equal 2, t[0].start
    assert_equal 7, t[0].stop
    assert_equal true, t[0].fragment?
    assert_equal 8, t[1].start
    assert_equal 16, t[1].stop
    assert_equal false, t[1].fragment?
    assert_equal 17, t[2].start
    assert_equal 22, t[2].stop
    assert_equal true, t[2].fragment?
  end
  
  
  def test_orf
    o = Orf::Orf.new
      
    o.aa_sequence = 'M*'
    assert_equal false, o.end_fragment?
      
    o.aa_sequence = 'KKK*'
    assert_equal true, o.end_fragment?
      
    o.aa_sequence = 'KKKT'
    assert_equal false, o.end_fragment?
  end
    
  # Am expecting to fail this for the moment, until I fix regex problems
  def test_middle_orf
    threads = OrfFinder.new.generate_longest_orfs('ATGTAGATGAAATAG')
    assert_equal 3, threads.length
    assert_equal 2, threads[0].length
    t = threads[0][0]
    assert_equal 0, t.start
    assert_equal 5, t.stop
    assert_equal false, t.fragment?
    t = threads[0][1]
    assert_equal 6, t.start
    assert_equal 14, t.stop
    assert_equal false, t.fragment?
  end
  
  
  def test_longest_orf
    finder = OrfFinder.new
    # test no orfs
    assert_nil finder.longest_orf('ATTTTTTT')
    
    # test 2 orfs in the same frame
    o = finder.longest_orf('ATGTAGATGAAATAG')
    assert o
    assert_equal 3, o.length
    
    # test 2 frames in different frames
    o = finder.longest_orf('ATGATGATGTAGAAAATGAAATAG')
    assert o
    assert_equal 4, o.length, o.aa_sequence
  end
  
  
  def test_longest_full_orf
    finder = OrfFinder.new
    # test no orfs
    assert_nil finder.longest_full_orf('ATTTTTTT')
    
    # test 2 orfs in the same frame, where the first is full
    o = finder.longest_full_orf('ATGTAGATGAAAAAAAAA')
    assert o
    assert_equal 'M*', o.aa_sequence
    assert_equal 2, o.length

    
    # test 2 frames in different frames
    o = finder.longest_full_orf('ATGATGATGTAGAAAATGAAATAG')
    assert o
    assert_equal 4, o.length, o.aa_sequence
  end
  
  def test_longest_m_orf
    finder = OrfFinder.new
    # test no orfs
    assert_nil finder.longest_m_orf('ATTTTTTT')
    
    # test 2 orfs in the same frame, where the first is full
    o = finder.longest_m_orf('ATGTAGATGAAAAAAAAA')
    assert o
    assert_equal 'MKKK', o.aa_sequence
    assert_equal 4, o.length

    
    # test 2 frames in different frames
    o = finder.longest_m_orf('ATGTAGAAAATGAAAAAA')
    assert o
    assert_equal 3, o.length, o.aa_sequence
  end
  
  def test_protein
    finder = OrfFinder.new
    
    o = finder.longest_protein_m_orf('M*')
    assert_equal 'M*', o.aa_sequence
    
    o = finder.longest_protein_m_orf('M*MEG')
    assert_equal 'MEG', o.aa_sequence
  end
  
  def test_command_line_protein
    proteins = system("orf_finder.rb -pf data/orf_finder_proteins.fa")
    expected = <<-END_OF_FASTA
>Contig_Oyster_90_8089_1_95_7630_1_98_6720_1
MLRFIAIVALIATVNAKGGTYGIGVLPSVTYVSGGGGGGGYYPGSYGTYGGGYPVTYGGFGPGSVYGSINSFGGVSTSAYGLYGTSPAVRGAAQGAATLSALGVASGVPSRVSGSSIGIGGGRALVSGSATPIGYYGVPYGGYSYGVPSYGYGYGYGYPSYGISYGYPGYGYGGYGGYGYPDVAHFGGSTYGNLATGAISSPTSGITIPYGGAIGLYGGYGGYGLGYGGYGGGYGGYGGGYGGYGGGYGGYYGGYYPSYGSSLTGVSQSLSFGRAVMSGQAFGAGVPAFGSVNFGNFGVGSGGIYGPSIYGGGGDIYGGGATIIX
>Contig_Oyster_90_3067_2_95_2874_2_98_2510_2
MLRFLAVVALIATVNAGGYGLYGGGGYPGIYGTYGGYPSIYGGFGPGGVYGSINSYGGVSTGAYGLYGTSPAVRGAAQGAATLSALGVASGVPSRVSGSSIGIGGGRALVSGSATPIGYYGVPYGGYSYGVPSYGYGYGYGYGYPSYGISYGYPGYGYGGYGGYGYPDVAYFGGSTYGNLATGAISSPTSGVTIPYGGALGLYGGYGGYGGYGLGYGGYGLGYGGYGLGYGGYGGGYGGYYGGYYPSYGSSLTGVSQSLSFGRAVMGGQAFGAGVPAFGSVNFGNFGVGTGGIYGPGIYGGGIYGGGGIYGGGATIIRRKKY*
    END_OF_FASTA
    assert_equal expected, proteins
  end
end