require "test/unit"
require 'tempfile'

class BlastEvalueGrepTest < Test::Unit::TestCase
  def test_default
    input = <<EOF
Aqu1.200001	Aqu1.200001	100.00	294	0	0	1	294	1	294	2e-166	 583
Aqu1.200001	Aqu1.200002	100.00	74	0	0	221	294	237	164	5e-35	 147
Aqu1.200001	Aqu1.200002	98.00	74	0	223	230	294	237	164	5e-4	 147
Aqu1.200002	Aqu1.200002	100.00	237	0	0	1	237	1	237	2e-132	 470 
EOF
    expected = <<EOF
Aqu1.200001	Aqu1.200001	100.00	294	0	0	1	294	1	294	2e-166	 583
Aqu1.200001	Aqu1.200002	100.00	74	0	0	221	294	237	164	5e-35	 147
Aqu1.200002	Aqu1.200002	100.00	237	0	0	1	237	1	237	2e-132	 470 
EOF
    Tempfile.open('test') do |t|
      t.puts input
      t.close
      output = `blast_evalue_grep.rb #{t.path}`
      assert_equal expected, output
    end
  end
  
  def test_none
    input = <<EOF
Aqu1.200001	Aqu1.200001	100.00	294	0	0	1	294	1	294	2e-166	 583
Aqu1.200001	Aqu1.200002	100.00	74	0	0	221	294	237	164	5e-35	 147
Aqu1.200001	Aqu1.200002	98.00	74	0	223	230	294	237	164	5e-34	 147
Aqu1.200002	Aqu1.200002	100.00	237	0	0	1	237	1	237	2e-132	 470 
EOF
    expected = <<EOF
Aqu1.200001	Aqu1.200001	100.00	294	0	0	1	294	1	294	2e-166	 583
Aqu1.200001	Aqu1.200002	100.00	74	0	0	221	294	237	164	5e-35	 147
Aqu1.200001	Aqu1.200002	98.00	74	0	223	230	294	237	164	5e-34	 147
Aqu1.200002	Aqu1.200002	100.00	237	0	0	1	237	1	237	2e-132	 470 
EOF
    Tempfile.open('test') do |t|
      t.puts input
      t.close
      output = `blast_evalue_grep.rb <#{t.path}`
      assert_equal expected, output
    end
  end
  
  def test_specified_evalue
    input = <<EOF
Aqu1.200001	Aqu1.200001	100.00	294	0	0	1	294	1	294	2e-166	 583
Aqu1.200001	Aqu1.200002	100.00	74	0	0	221	294	237	164	5e-35	 147
Aqu1.200001	Aqu1.200002	100.00	74	0	0	221	294	237	164	5e-34	 147
Aqu1.200001	Aqu1.200002	98.00	74	0	223	230	294	237	164	5e-4	 147
Aqu1.200002	Aqu1.200002	100.00	237	0	0	1	237	1	237	2e-132	 470 
EOF
    expected = <<EOF
Aqu1.200001	Aqu1.200001	100.00	294	0	0	1	294	1	294	2e-166	 583
Aqu1.200001	Aqu1.200002	100.00	74	0	0	221	294	237	164	5e-35	 147
Aqu1.200002	Aqu1.200002	100.00	237	0	0	1	237	1	237	2e-132	 470 
EOF
    Tempfile.open('test') do |t|
      t.puts input
      t.close
      output = `blast_evalue_grep.rb -e 5e-35 #{t.path}`
      assert_equal expected, output
    end
  end
end
