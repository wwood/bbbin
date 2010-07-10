require "test/unit"
require 'tempfile'
# These examples are sometimes taken from http://genome.jgi-psf.org/Triad1/Triad1.download.html
class BlastSaturationTest < Test::Unit::TestCase
  # create blast databases
  def create_formatted_blast_database(fasta_string)
    @db = Tempfile.new('db')
    @db.puts fasta_string
    @db.close
    `formatdb -i #{@db.path}`
  end
  
  def test_no_hits
    create_formatted_blast_database <<EOF
>seq1
ATGCATGCATGC
>seq2
AGAGAGAGAGAGAGA
>jgi|Triad1|62822|fgeneshTA2_pg.C_scaffold_1150000001_modified
MSDIKKFKKELDTFHNSEATRKSSNKIINVICAKSDNYIGGSADLSGSNGTKCN
EOF

    Tempfile.open('fasta') do |t|
      t.puts '>notHitting'
      t.puts 'GGGGGGGGG'
      t.close
      Tempfile.open('output') do |out|
        `blast_saturation.rb #{t.path} '-p blastp -d #{@db.path} -F F -z 10000' >#{out.path}`
        output = `cat #{out.path}`
        assert_equal '', output
      end
    end
  end
  
  def test_one_hit
    create_formatted_blast_database <<EOF
>seq1
ATGCATGCATGC
>seq2
AGAGAGAGAGAGAGA
>jgi|Triad1|62822|fgeneshTA2_pg.C_scaffold_1150000001_modified
MSDIKKFKKELDTFHNSEATRKSSNKIINVICAKSDNYIGGSADLSGSNGTKCN
EOF
    
    Tempfile.open('fasta') do |t|
      t.puts '>jgi|Triad1|62822|fgeneshTA2_pg.C_scaffold_1150000001'
      t.puts 'MSDIKKFKKELDTFHNSEATRKSSNKIINVICAKSDNYIGGSADLSGSNGTKCNNHKIISKEDFTGNYIHYGVREHAMVA'
      t.close
      Tempfile.open('output') do |out|
        # use -z 100 so that the db can be changed without minor e-value changes
        `blast_saturation.rb #{t.path} '-p blastp -d #{@db.path} -F F -z 100' >#{out.path}`
        output = File.open(out.path).read
        assert_equal [
        'jgi|Triad1|62822|fgeneshTA2_pg.C_scaffold_1150000001',
        'jgi|Triad1|62822|fgeneshTA2_pg.C_scaffold_1150000001_modified',
        '100.00',
        '54',
        '0',
        '0',
        '1',
        '54',
        '1', '54',  '9e-30',  ' 109'].join("\t"), output.strip
      end
    end
  end
  
  
  def test_three_hits
    create_formatted_blast_database <<EOF
>seq1
ATGCATGCATGC
>seq2
AGAGAGAGAGAGAGA
>jgi|Triad1|62822|fgeneshTA2_pg.C_scaffold_1150000001_modified
MSDIKKFKKELDTFHNSEATRKSSNKIINVICAKSDNYIGGSADLSGSNGTKCN
>third
KSFGGTFLVFSD
>second
FFLVSLSLASIPAVNMAMYLEKKVSYLYEASDGF
EOF
    
    Tempfile.open('fasta') do |t|
      t.puts '>jgi|Triad1|62822|fgeneshTA2_pg.C_scaffold_1150000001'
      t.puts 'MSDIKKFKKELDTFHNSEATRKSSNKIINVICAKSDNYIGGSADLSGSNGTKCNNHKIISKEDFTGNYIHYGVREHAMVAIMNGINIHDKFKSFGGTFLVFSDFFLVSLSLASIPAVNMAMYLEKKVSYLYEASDGF'
      t.close
      
      output = `blast_saturation.rb #{t.path} '-p blastp -d #{@db.path} -F F -z 10000 -b 1 -v 1'`
      expected = <<EOF
jgi|Triad1|62822|fgeneshTA2_pg.C_scaffold_1150000001	jgi|Triad1|62822|fgeneshTA2_pg.C_scaffold_1150000001_modified	100.00	54	0	0	1	54	1	54	6e-28	 110
jgi|Triad1|62822|fgeneshTA2_pg.C_scaffold_1150000001	second	100.00	34	0	0	104	137	1	34	2e-15	68.9
jgi|Triad1|62822|fgeneshTA2_pg.C_scaffold_1150000001	third	100.00	12	0	0	92	103	1	12	0.007	26.9
EOF
      assert_equal expected, output
    end
  end
end
