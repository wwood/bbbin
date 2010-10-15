require "test/unit"
require 'mask_eupathdb_domains'

class EuPathDBDomainAnnotationFileTest < Test::Unit::TestCase
  def test_parse_line
    line = <<EOF
PFI1830c	SUPERFAMILY	SSF101967	0042205	1895	2072	  1.7E-03
EOF
    obs = EuPathDB::DomainAnnotationFile.new(line)
    domains = obs.domain_hash
    assert_kind_of Hash, domains
    assert domains['PFI1830c']
    assert_kind_of Array, domains['PFI1830c']
    assert_equal 1, domains['PFI1830c'].length
    d = domains['PFI1830c'][0]
    assert_equal 'PFI1830c', d.protein_accession
    assert_equal 1895, d.start
    assert_equal 2072, d.stop
    assert_equal 1.7e-3, d.evalue
  end
  
  def test_parse_file
    obs = EuPathDB::DomainAnnotationFile.new(File.open(File.join('data','plasmodb_domain_definition_extract.txt')))
    domains = obs.domain_hash
    assert_equal 2, domains.keys.length
    assert_equal 2, domains['PFI1830c'].length
    assert_equal 1, domains['PFI1825w'].length
    assert_equal 1.7E-03, domains['PFI1830c'][1].evalue
    assert_equal 0.0, domains['PFI1825w'][0].evalue
  end
end
