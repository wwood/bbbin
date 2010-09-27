require "test/unit"
require 'tempfile'
# These examples are sometimes taken from http://genome.jgi-psf.org/Triad1/Triad1.download.html
class BlastSaturationTest < Test::Unit::TestCase
  #create blast databases
  def create_formatted_blast_database(fasta_string)
    @db = Tempfile.new('db')
    @db.puts fasta_string
    @db.close
    `formatdb -i #{@db.path}`
  end
  
  #create blast+ databases
  def create_formatted_blast_plus_database(fasta_string)
    @db = Tempfile.new('db')
    @db.puts fasta_string
    @db.close
    `makeblastdb -in #{@db.path}`
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
        `blast_saturation.rb -q -f #{t.path} -b '-p blastp -d #{@db.path} -F F -z 10000' >#{out.path}`
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
        #use -z 100 so that the db can be changed without minor e-value changes
        `blast_saturation.rb -q -f #{t.path} -b '-p blastp -d #{@db.path} -F F -z 100' >#{out.path}`
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
      
      output = `blast_saturation.rb -q -f #{t.path} -b '-p blastp -d #{@db.path} -F F -z 10000 -b 1 -v 1'`
      expected = <<EOF
jgi|Triad1|62822|fgeneshTA2_pg.C_scaffold_1150000001	jgi|Triad1|62822|fgeneshTA2_pg.C_scaffold_1150000001_modified	100.00	54	0	0	1	54	1	54	6e-28	 110
jgi|Triad1|62822|fgeneshTA2_pg.C_scaffold_1150000001	second	100.00	34	0	0	104	137	1	34	2e-15	68.9
jgi|Triad1|62822|fgeneshTA2_pg.C_scaffold_1150000001	third	100.00	12	0	0	92	103	1	12	0.007	26.9
EOF
      assert_equal expected, output
    end
  end
  
  def test_blast_plus_bad_program
    Tempfile.open('err') do |t|
      `blast_saturation.rb -p notaBLASTprogram '-d nowhere' 2>#{t.path}`
      expected = <<EOF
Unexpected blast+ program found: `notaBLASTprogram'. Expected one of blastn, blastp, blastx, tblastn, tblastx.
EOF
      assert_equal expected, t.read
    end
  end
  
  def test_blast_plus_three_hits
    create_formatted_blast_plus_database <<EOF
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
      
      output = `blast_saturation.rb -q -p blastp -f #{t.path} -b '-db #{@db.path} -seg no -dbsize 10000 -num_descriptions 1 -num_alignments 1'`
      expected = <<EOF
jgi|Triad1|62822|fgeneshTA2_pg.C_scaffold_1150000001	jgi|Triad1|62822|fgeneshTA2_pg.C_scaffold_1150000001_modified	100.00	54	0	0	1	54	1	54	6e-28	 110
jgi|Triad1|62822|fgeneshTA2_pg.C_scaffold_1150000001	second	100.00	34	0	0	104	137	1	34	2e-15	68.9
jgi|Triad1|62822|fgeneshTA2_pg.C_scaffold_1150000001	third	100.00	12	0	0	92	103	1	12	0.007	26.9
EOF
      assert_equal expected, output
    end    
  end
  
  # Blastx is different because the masking was previously failing. Actually,
  # turns out the reverse bit was the problem.
  def test_blastx_reverse
    create_formatted_blast_plus_database <<EOF
>gi|164421170|ref|YP_001648654.1| cytochrome c oxidase subunit III [Hippospongia lachne]
MKYYRPYHLVDSSPWPFLGGCAGLSLVVGGILYMHYGYIWLMVSGVLFVGIIMVVWWRDVVRESTFLGKH
NTIVKRGIKYGMILFIVSEIMLFFSLFWAFFHNSLSFAVELGGCWVPRGVEALDYKAVPLLNTALLLGSG
MLVTWAHYGIIRGWRVVGLRALGSGVALGLLFTCLQGFEYYVASFTIADSVYGSVFYLMTGAHGLHVIIG
SVFLTVCWFRLLYYQFRDDDPVGFELAAWYWHFVDVVWLFLYIFVYCWGS
EOF
    Tempfile.open('fasta') do |t|
      t.puts '>Carteriospongia_foliascens_Contig20_chopped_2000-2090'
      t.puts 'CCCCAAAAAGGGCCATGGGCTAGAGTCTACTAAATGATACGACCTATAATATTTCATTCCCCAAACATAAAAAATATTGATATTATTTTGT'
      t.close
      
      output =  `blast_saturation.rb -q -p blastx -f #{t.path} -b '-db #{@db.path} -seg no -dbsize 10000 -num_descriptions 1 -num_alignments 1'`
      expected = <<EOF
Carteriospongia_foliascens_Contig20_chopped_2000-2090	gi|164421170|ref|YP_001648654.1|	94.74	19	1	0	57	1	1	19	6e-09	45.4
EOF
      assert_equal expected, output
    end 
  end
end
