require "test/unit"
require 'mask_eupathdb_domains'

class EuPathDBDomainAnnotationFileTest < Test::Unit::TestCase
  
  def test_parse_line
    line = <<EOF
PFI1830c  PFAM  PF05658 Hep_Hag 1894  1921    3.3E-05
EOF
    obs = EuPathDB::DomainAnnotationFile.new
    assert_kind_of EuPathDB::Domain, obs
    assert_equal 'PFI1830c', obs.protein_accession
    assert_equal 1894, obssstart
    assert_equal 1921, obs.stop
  end
end